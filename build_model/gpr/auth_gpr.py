import logging
import os
import re
from io import StringIO
from typing import List, Optional, Tuple

import pandas as pd
import requests

from gpr.ast_gpr import sanitize_gpr

LOGGER = logging.getLogger(__name__)
GPRURL2 = "http://www.genome.jp/dbget-bin/www_bget?ec:"
REMOVE_PAT = re.compile(r"\)|\(|'| |\]|\[")
TRANSFERRED_PAT = re.compile(
    r'It is now listed as..\n+.+EC\-([0-9\.\-]+)" class\="EC\-NUMBER"', re.MULTILINE
)


def read_env_biocyc():
    """Get the email and password from `BIOCYC_*` env variables.

    BIOCYC_EMAIL: email of login user
    BIOCYC_PASSWORD: the password

    To have this, you must create an account in BioCyc, given that you are
    connecting from a licensed institution.
    """

    email = os.environ["BIOCYC_EMAIL"]
    password = os.environ["BIOCYC_PASSWORD"]
    return email, password


def setup_biocyc_session(
    email: Optional[str] = None, password: Optional[str] = None
) -> requests.Session:
    s = requests.Session()  # create session
    # Post login credentials to session:
    if email is None or password is None:
        email, password = read_env_biocyc()
    s.post(
        "https://websvc.biocyc.org/credentials/login/",
        data={"email": email, "password": password},
    )
    return s


def get_ecnumber_biocyc_html(
    ec_number: str, session: requests.Session, org: str = "META"
) -> str:
    """Get raw HTML representing a EC entry in BioCyc."""
    return session.get(
        f"https://websvc.biocyc.org/{org}/NEW-IMAGE?type=EC-NUMBER&object=EC-{ec_number}"
    ).text


def get_html(request_url: str, session: Optional[requests.Session] = None) -> str:
    """Fetch an html by perfoming a GET HTTPS request, maybe with session."""
    if session is not None:
        return session.get(request_url).text
    else:
        return requests.get(request_url).text


def pattern_match_org(page: str, org: str = "Homo Sapiens") -> List[str]:
    return sorted(
        [
            x[0]
            for x in [
                re.findall(
                    r"\<b\>Gene:\</b\>..([a-zA-Z0-9]+)....([A-Za-z0-9:,-_\/]+)<br>.*"
                    + org,
                    str(x),
                )
                for x in page.splitlines()
            ]
            if x
        ]
    )


def match_biocyc_page(page: str, humancyc: bool = False) -> List[str]:
    if humancyc:
        # since HumanCyc relates to human,
        # we don't need to look for a particular organism
        return pattern_match_org(page, "")
    # fetch in human
    urls0 = pattern_match_org(page, "Homo Sapiens")
    # if failed, failed in mouse
    if not urls0:
        urls0 = pattern_match_org(page, "Mus Musculus")
    return urls0


def getGPR(
    ec_number: str, session: Optional[requests.Session]
) -> Optional[Tuple[List[str], str, str, str, str]]:
    """Retrieve GPR given EC-number for BioCyc."""
    if session is None:
        session = setup_biocyc_session()
    page = get_ecnumber_biocyc_html(ec_number, session, org="HUMAN")
    LOGGER.info(f"Page length: {len( page )}")
    urls0 = match_biocyc_page(page, humancyc=True)
    LOGGER.info(f"Urls matched from matched from HUMAN: {len(urls0)}")
    # if failed, maybe the EC number was changed
    if not urls0:
        ec_number_match = TRANSFERRED_PAT.findall(page)
        if ec_number_match:
            LOGGER.info(f"Trying transferred EC number: {ec_number_match[0]}")
            gprs = getGPR(ec_number_match[0], session)
            if gprs is not None:
                return gprs
        else:
            pass
    # if failed, try MetaCyc
    if not urls0:
        page = get_ecnumber_biocyc_html(ec_number, session)
        urls0 = match_biocyc_page(page)
        LOGGER.info(f"Urls matched from changed MetaCyc: {len(urls0)}")
    if not urls0:
        LOGGER.warn("Trying Kegg")
        # return fetch_kegg_rest(ec_number)
        parsed = _fetch_kegg_from_ec_html(ec_number)
    else:
        parsed = parseGPR(urls0, page, session)
    if parsed[0]:
        a, b, c, d, gpr = parsed
        gpr = sanitize_gpr(gpr)
        gpr = "([" + re.sub(r"([A-Z0-9\-\.]+)", r"([\1])", gpr) + "])"
        parsed = a, b, c, d, gpr 
    return parsed


def _fetch_kegg_from_ec_html(ec_number: str):
    #GPRPage2 = get_html(str(GPRURL2) + str(ec_number))
    GPRPage2 = get_html(str(GPRURL2) + str(ec_number))
    LOGGER.info(f"Page length from kegg: {len(GPRPage2)}")
    urls0 = sorted(
        set(
            REMOVE_PAT.sub(
                "",
                str(
                    re.findall(
                        r"(\([A-Za-z0-9]+\))",
                        str(
                            re.findall(
                                r"hsa:............................................",
                                GPRPage2,
                            )
                        ),
                    )
                ),
            ).split(",")
        )
    )
    LOGGER.info(f"Urls from HUMAN kegg: {urls0}")
    if not urls0[0]:  # if not for homo sapiens try with mus musculus
        urls0 = sorted(
            set(
                REMOVE_PAT.sub(
                    "",
                    str(
                        re.findall(
                            r"(\([A-Za-z0-9]+\))",
                            str(
                                re.findall(
                                    r"mmu:............................................",
                                    GPRPage2,
                                )
                            ),
                        )
                    ),
                ).split(",")
            )
        )
        LOGGER.warn(f"Urls from MOUSE kegg: {urls0}")
    if urls0[0]:
        urls1 = [x[0:] for x in urls0]
        urls2 = urls1
        urls3 = (
            "[(["
            + str(urls1)
            .replace("[", "")
            .replace("]", "")
            .replace(", ", "*1]) or ([")
            .replace("'", "")
            + "*1])]"
        )
        urls4 = re.sub(r"*[0-9]+", "", str(urls3))
    else:
        urls1 = ""
        urls2 = ""
        urls3 = ""
        urls4 = ""
    if not urls1:
        return ([], '', '', '', '')
    return urls0, urls1, urls2, urls3, urls4


def fetch_kegg_rest(ec_number: str) -> Optional[Tuple[List[str], str, str, str, str]]:
    """Fetch ec-number from kegg and link it to its genes.

    Perform the following requests:

    .. code-block::

        LINK ec-number -> KO (Kegg orthology)
        LINK KO (Kegg orthology) -> genes

    We map first to KO because LINK to genes only works in this way.

    There is also the possibility of linking from ec-number to HSA
    (human gene Kegg ids) the problem is that gene -> gene symbol
    requires additional parsing.
    """
    ec_to_ko = pd.read_csv(
        StringIO(get_html(f"https://rest.kegg.jp/link/ko/ec:{ec_number}")),
        sep="\t",
        names=["ec", "ko"],
    )
    LOGGER.info(f"{len(ec_to_ko)} Kegg orthologies fetched.")

    ko_to_genes = pd.read_csv(
        StringIO(
            get_html(
                f"https://rest.kegg.jp/link/genes/{'+'.join(ec_to_ko.ko.to_list())}"
            )
        ),
        sep="\t",
        names=["ko", "genes"],
    )

    genes = ko_to_genes.loc[
        ko_to_genes.genes.str.startswith("hsa:") & ~ko_to_genes.genes.isnull(), "genes"
    ].to_list()
    if not genes:
        return None
    LOGGER.info(f"{len(genes)} Kegg genes fetched.")
    genes_stripped = [gene.replace("hsa:", "") for gene in genes]
    return (
        genes,
        genes,
        genes_stripped,
        "[" + " or ".join(f"({genes_stripped}) *1") + "]",
        f"[{' or '.join(genes_stripped)}]",
    )


def parseGPR(
    urls0: List[str], page: str, session: requests.Session
) -> Tuple[List[str], str, str, str, str]:
    urls1 = [
        i for i in reversed(sorted([x[0:][0] for x in urls0], key=len))
    ]  # Sort genes by name lenght
    Ls = list("ABCDEFGHIJKLMNOPQRSTUVWXYZ")
    urls2 = []
    # urls22 = []
    for x in urls1:
        for y in urls0:
            if re.findall(r"\('" + str(x) + "',", str(y)):
                urls2 = urls2 + [y[1]]
                # urls22 = urls22 + [(y[0], y[1], Ls.pop(0))]
    # LOGGER.debug(f"{urls22=}")
    a = page.replace("\n", "").replace("Enzymes and Genes:", "\nEnzymes and Genes:")
    a2 = re.findall(r"Enzymes and Genes:.*[\S\s]+", a)
    b = a2[0].replace("<br> <a href=", "\n<br> <a href=").replace("</a>)", "</a>)\n")
    c = re.findall(
        r"<br> <a href=.*[\S\s].*>Homo sapiens</i>", b
    )  # evaluate the ezymes active in human Extracted from Biocyc# c = re.findall(r'<br> <a href=.*[\S\s].*>Homo sapiens</a>',b)
    if not c:
        c = re.findall(
            r'href="/gene\?orgid=HUMAN&id=(.*?)</a>', b
        )  # evaluate the ezymes active in human extracted from HumanCyc# c = re.findall(r'<br> <a href=.*[\S\s].*>Mus musculus</a>',b)
    if not c:
        c = re.findall(
            r"<br> <a href=.*[\S\s].*>Mus musculus</i>", b
        )  # evaluate the ezymes active in human # c = re.findall(r'<br> <a href=.*[\S\s].*>Mus musculus</a>',b)
    dict = {
        "<SUB>": "*",
        "<@SUB>(": " and (",
        "<@SUB>[": " and [",
        "][": "] and [",
        ")(": " and ",
        "'": "",
        ",": "",
        ";": "",
        "<BR>": "",
    }
    dict2 = {
        ", ": "*",
        "), (": " and ",
        "(": "",
        ")": "",
        ",": "",
        ";": "",
        "'": "",
        "<b>": "",
        "</b>": "",
        "<br>": "",
        "</br>": "",
        "<i>": "",
        "</i>": "",
        "<a>": "",
        "</a>": "",
        ":": "",
        "=": ",",
    }
    gpr2 = []
    i = 0
    while i < len(c):  # isoforms
        if re.findall(r'class="GENE"', c[i]):
            c2 = c[i].replace(" ", "\n")
            c3 = c2.replace('"\n', '" \n')
            isourl = re.findall(
                r"(http://biocyc.org/META/NEW-IMAGE\?type=ENZYME.*[\S\s])\" \n",
                c3.replace("/META", "http://biocyc.org/META"),
            )  # if it is a complex
            if not isourl:
                isourl = re.findall(
                    r"<br>\n<a\nhref=\"(http://biocyc.org/gene\?orgid.*)\"",
                    c3.replace("/gene?orgid", "http://biocyc.org/gene?orgid"),
                )  # if it is  not a complex
            # if not isourl: isourl = ['https://biocyc.org/gene?orgid=META&id='+re.findall('^(.*?)"',c3)[0].strip()] # if it is taken from HumanCyc
            if not isourl:
                isourl = [
                    "https://websvc.biocyc.org/gene?orgid=HUMAN&id="
                    + re.findall('^(.*?)"', c3)[0].strip()
                ]  # if it is taken from HumanCyc
            isopage = str(
                get_html(isourl[0].replace(" ", "").replace('"', ""), session)
            )
            d1 = isopage.replace("\n", " ").replace(
                "</a>", "\n</a>"
            )  # d1 =  isopage.replace('\n',' ').replace('<br>','<br>\n')
            d2 = [
                x.strip() for x in re.findall(r"Subunit composition.*=(.*)\\n", d1)
            ]  # re.findall(r'</a> </p>  <p class=ecoparagraph>  Subunit composition .*=(.*)', d1.replace(" <br>",""))
            if not d2:
                d2 = re.findall(
                    r"</a> </p>  <p class=ecoparagraph>  Subunit Composition.*=(.*)",
                    d1.replace(" <br>", ""),
                )
            if not d2:
                d2 = re.findall(r"Subunit composition.*=(.*)", d1.replace(" <br>", ""))
                if d2:
                    d2 = d2[0]
            # if not d2:
            #    d2 = re.findall(r'Subunit Composition.*=(.*)', d1.replace(" <br>",""))
            #    if d2: d2 = d2[0]
            if d2:
                d2 = d2[
                    0
                ]  # if we found subunit composition then only use the first match
            if not d2:
                isourl2 = re.findall(r"/gene-tab.*META&tab=SUMMARY", str(isopage))
                if isourl2:
                    isourl2 = ["https://websvc.biocyc.org" + isourl2[0]]
                    isopage2 = str(get_html(isourl2[0], session))
                    d2 = re.findall(
                        r"Subunit Composition</td><td align=LEFT valign=TOP>(.*)",
                        isopage2.replace("<tr><td", "\n").replace("</td></tr>", "\n"),
                    )
            if not d2:
                d1 = str(isopage).replace("\n", " ").replace("</a>", "</a>\n")
                d2 = re.findall(
                    "/META/NEW-IMAGE\?type=REACTION.*\</a>",
                    d1.replace(" <br>", ""),
                )
            if not d2:
                # isourl = 'https://biocyc.org/gene?orgid=META&id='+re.findall('^(.*?)"',c3)[0].strip() # if it is taken from HumanCyc
                isourl = (
                    "https://websvc.biocyc.org/gene?orgid=HUMAN&id="
                    + re.findall('^(.*?)"', c3)[0].strip()
                )  # if it is taken from HumanCyc
                isopage = str(get_html(isourl, session))
                # isopage = getHtml(isourl)
                # used_ip = []
                # isopage = getHtml3(url, used_ip, 10)
                # isopage = getHtml_non_anonymous_2(isourl)
                # isopage = getHtml_anonymous_2(isourl, used_ip, memory_limit)
                d1 = str(isopage).replace("\n", " ").replace("</a>", "\n</a>")
                d2 = re.findall(
                    r">Subunit Composition</td><td.*?>(.*?)</td>", str(isopage)
                )
                if not d2:
                    ith_gene = re.findall(
                        r"gene\\n\\n<br><font class=\\\'header\\\'>(.*?)<",
                        isopage,
                    )
                    if ith_gene:
                        ith_gene = ith_gene[0].strip()
                        position = str(isopage).find("[" + ith_gene + "]<SUB>")
                        if position == -1:
                            d2 = ith_gene
                        else:
                            d2 = re.findall(
                                r"^(.*?</SUB>)",
                                isopage[position : position + 100],
                            )[0].strip()
            LOGGER.debug(f"{d2=}")
            # if re.findall(r'<br> <b>Synonyms:</b>.*><b>(.*)<', d1): synonim = [multiple_replace(dict2, multiple_replace(dict,x)) for x in re.findall(r'<br> <b>Synonyms:</b>.*><b>(.*)<', d1)]
            synonim = [
                multiple_replace(dict2, multiple_replace(dict, x))
                for x in re.findall(r"<br> <b>Synonyms:</b>.*><b>(.*)<", d1)
            ]
            if not synonim and re.findall(r'top">\\nSynonyms</td>', d1):
                synonim = [
                    multiple_replace(dict2, multiple_replace(dict, x))
                    for x in re.findall(r'Synonyms</td>.*top">(.*?)</td></tr>', d1)
                ]  # in case it comes from HumanCyc using getHtml2 function
                if synonim:
                    synonim = synonim[0].strip().split(";")
            if synonim:
                synonim = list(
                    filter(
                        None,
                        [re.sub(" $", "", re.sub("^ ", "", x)) for x in synonim],
                    )
                )
            d = list(
                filter(
                    None,
                    re.findall(
                        "[\S]+",
                        str(d2)
                        .upper()
                        .replace("</SUB>", "<@SUB>")
                        .replace("/", " ")
                        .replace("(", "[")
                        .replace(")", "]"),
                    ),
                )
            )  # original without "list"
            j = 0
            e = [""] * len(d)
            while j < len(d):
                z = 0
                while z < len(urls1):
                    IsVar = ""
                    if not re.findall("/NEW-IMAGE\?TYPE=REACTION", d[j].upper()):
                        if re.findall(urls1[z].upper(), d[j].upper()):
                            IsVar = [d[j].upper()]
                        if not IsVar and synonim:
                            dict3 = {}
                            for x in synonim:
                                dict3[x.upper()] = urls1[z].upper()
                            IsVar = list(
                                filter(
                                    None,
                                    [
                                        re.findall(
                                            x.upper()
                                            .replace("(", "\(")
                                            .replace(")", "\)")
                                            .replace('"', '"')
                                            .replace("'", "'"),
                                            d[j].upper(),
                                        )
                                        for x in synonim
                                    ],
                                )
                            )  # original without list
                            if IsVar:
                                IsVar = [multiple_replace(dict3, x) for x in IsVar[0]]
                            del dict3
                        if not IsVar:
                            IsVar = re.findall(urls1[z].upper(), d[j].upper())
                    if IsVar:
                        h = (
                            multiple_replace(dict, IsVar[0])
                            .replace("<@SUB>", "")
                            .replace("] and [", "]*1 and [")
                            .replace("] or [", "]*1 or [")
                            .replace(")", "")
                            .replace("(", "")
                        )
                        h = re.sub("]$", "]*1", h)
                        if not re.findall(r"]", h) or not re.findall(r"\[", h):
                            h = (
                                "[" + str(h) + "]*1"
                            )  # to adapt the case of single gene reaction association
                        h = "[" + str(h) + "]*1"
                        if re.findall(r"Subunit composition.*=", d1) or re.findall(
                            r"Subunit Composition:.*=", d1
                        ):
                            z = len(
                                urls1
                            )  # if we are analyzing a complex, then if the gpr is found stop the iteration
                            s = h
                        else:
                            NestLevel = re.findall(r"(?=(\]\*[0-9]+\]\*))", h)
                            h0 = ParseNestedParen(str(h), len(NestLevel))[0]
                            Mult = re.findall(r"([\]\*0-9]+$)", h)
                            if Mult:
                                Mult = eval(
                                    str(Mult[0].split("*")[1:])
                                    .replace("'", "")
                                    .replace("]", "")
                                    .replace("[", "")
                                    .replace(", ", "*")
                                )
                            else:
                                Mult = 1
                            h = h0
                            for x in urls1:
                                h = h.replace(x, "@")
                            h = filter(None, h.split("@"))
                            dict4 = {
                                "\\": "",
                                '"': "",
                                "'": "",
                                "?": "",
                                ")": "",
                                "(": "",
                                "]": "",
                                "[": "",
                                "*": "",
                                "$": "",
                                ".": "",
                            }
                            hh = []
                            for x in h:
                                for w in urls1:
                                    for y in urls1:
                                        x1 = multiple_replace(dict4, x)
                                        m1 = str(y) + x1 + str(w)
                                        m2 = str(y) + x1
                                        if re.findall(
                                            m1, multiple_replace(dict4, h0)
                                        ) or re.findall(
                                            m2 + "$",
                                            multiple_replace(dict4, h0),
                                        ):
                                            t = re.findall(r"([0-9]+)", str(x))
                                            if t:
                                                t = str(
                                                    eval(
                                                        str(t)
                                                        .replace("'", "")
                                                        .replace('"', "")
                                                        .replace("]", "")
                                                        .replace("[", "")
                                                        .replace(", ", "*")
                                                    )
                                                )
                                            else:
                                                t = "1"
                                            hh = hh + [str(y) + "*" + t]
                            hh = sorted(set(filter(None, hh)))
                            if not hh:
                                hh = [
                                    "[" + str(multiple_replace(dict4, h0)) + "]*1"
                                ]  # no hh means that the string only contains the gene name # if not hh: hh = [str(h0)+"*1"]
                            h22 = str(hh).replace("'", "").replace(", ", " and ")
                            UniqGene = hh  # This line has been added and the next has been commented to express in a SGPR dimers, trimers, tetramers, ....
                            # UniqGene = sorted(set([x.split("*")[0] for x in hh])) # if some of the subunits of the complex encoded by the same gene it is necessary to change"S and S" by 2*A # CHECK TO ADAPT TO THE GRASP CASE
                            GeneCount = [
                                (x.split("*")[0], x.split("*")[1])
                                for x in h22.replace("[", "")
                                .replace("]", "")
                                .split(" and ")
                            ]
                            UniqCount = []
                            for x in UniqGene:
                                count = "0"
                                for y in GeneCount:
                                    if re.findall(x, y[0]):
                                        count = count + "+" + y[1]
                                UniqCount = UniqCount + [[x, str(eval(count))]]
                            s = (
                                "["
                                + str([str(x).replace("', '", "*") for x in UniqCount])
                                .replace("'", "")
                                .replace('"', "")
                                .replace(", ", " and ")
                                + "]*"
                                + str(Mult)
                            )
                            z = len(urls1)
                    else:
                        s = ""
                        z += 1
                e[j] = s  # multiple_replace(dict2, s)
                j += 1
            e = sorted(set(filter(None, e)))  # remove empty and duplicate elements
            gpr = (
                str(e)
                .replace("'[", "")
                .replace("]'", "")
                .replace("'", "")
                .replace(", ", " and ")
            )
            gpr2 = gpr2 + [
                str(gpr).replace("'[", "").replace("]'", "").replace("'", "")
            ]
        i += 1
    urls3 = [
        str(sorted(set(gpr2)))
        .replace("'[", "[")
        .replace("]'", "]")
        .replace(", ", "] or [")
        .replace("'", "")
        .replace("*0", "*1")
    ]
    testIs = list(filter(None, [x.replace("[", "").replace("]", "") for x in urls3]))
    if not urls3 or not testIs:
        urls3 = [str(x) + "*1" for x in urls1]
    g = []
    r = 0
    while r < len(urls3):
        GsprGenes = sorted(
            set(
                re.findall(
                    r"[A-Za-z0-9]+",
                    re.sub(
                        r"\*[0-9]+",
                        "",
                        str(urls3[r]).replace(" and ", " ").replace(" or ", " "),
                    ),
                )
            )
        )
        GsprGenes = [x.upper().replace("[", "").replace("]", "") for x in GsprGenes]
        intersection = len(GsprGenes) - len(
            set(GsprGenes).intersection([x.upper() for x in urls1])
        )
        if intersection < 1:
            g = g + [
                str(urls3[r].upper()).replace(" AND ", " and ").replace(" OR ", " or ")
            ]
        r += 1
    urls3 = (
        str(g)
        .replace("'', ", "")
        .replace("', ''", ")")
        .replace("', '", ") or (")
        .replace("'", "(")
        .replace("](]", "])]")
        .replace('"', "")
        .replace(" *", "*")
        .replace("[ ", "[")
        .replace(" ]", "]")
        .replace("(]", ")]")
        .replace("[)", "[(")
    )
    if re.findall(r"Subunit composition.*=", d1):
        urls3 = urls3
    else:
        w = []
        for x in urls3.split(" or "):  # Evaluate Isoforms
            if re.findall(r"[A-Za-z0-9].*?", x):  # eliminate empty spots
                dict4 = {}
                for y in x.split(" and "):  # Evaluate complex
                    if re.findall(r"[A-Za-z0-9].*?", y):  # eliminate empty spots
                        y2 = re.sub("(\(|\)|\]|\[)", "", y)
                        a = re.findall(r"([A-Za-z0-9\-]+)\*([0-9\*]+)", y2)
                        LOGGER.debug(f"{a=}")
                        if not a[0][0] in dict4:
                            dict4[a[0][0]] = eval(a[0][1])
                        else:
                            dict4[a[0][0]] = eval(str(dict4[a[0][0]]) + "+" + a[0][1])
                        LOGGER.debug(f"{dict4=}")
                z = [(a + "*" + str(dict4[a])) for a in dict4]
                w = w + [z]
        urls3 = (
            str(
                sorted(set([str(x).replace(", ", " and ").replace("'", "") for x in w]))
            )
            .replace("'", "")
            .replace(", ", " or ")
            .replace("[", "([")
            .replace("]", "])")
        )
    urls4 = re.sub(r"\*[0-9]+", "", urls3)
    if not urls4: 
        return ([], '', '', '', '')
    else:
        return urls0, urls1, urls2, urls3, urls4


def multiple_replace(dict, text):
    """ "Multiple Replacement."""
    # Create a regular expression  from the dictionary keys
    regex = re.compile("(%s)" % "|".join(map(re.escape, dict.keys())))
    # For each match, look-up corresponding value in dictionary
    return regex.sub(lambda mo: dict[mo.string[mo.start() : mo.end()]], text)


def ParseNestedParen(string, level):
    """Generate strings contained in nested (), indexing i = level"""
    if len(re.findall(r"\[", string)) == len(re.findall(r"\]", string)):
        LeftRightIndex = [
            x
            for x in zip(
                [Left.start() + 1 for Left in re.finditer(r"\[", string)],
                reversed([Right.start() for Right in re.finditer(r"\]", string)]),
            )
        ]
    elif len(re.findall(r"\[", string)) > len(re.findall(r"\]", string)):
        return ParseNestedParen(string + "]", level)
    elif len(re.findall(r"\[", string)) < len(re.findall(r"\]", string)):
        return ParseNestedParen("[" + string, level)
    else:
        return "fail"
    return [string[LeftRightIndex[level][0] : LeftRightIndex[level][1]]]
