from functions_create_figure import *

# Load Models

model1 = cobra.io.read_sbml_model('models/Human-GEM_2022-06-21.xml')
model12 = cobra.io.read_sbml_model('models/THG-beta1.xml')
model15 = cobra.io.read_sbml_model('models/THG-beta2.xml')
modelTHG = cobra.io.read_sbml_model('models/THG-2023-02-25.xml')

file = 'reports/THG-beta-1_vs_Human1.xlsx'
# file = 'reports/THG-beta-2_vs_Human1.xlsx'
# file = 'reports/THG_vs_Human1.xlsx'


# Compare metabolite ids (model1 vs model12) 

sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})

df=pd.read_excel(file,sheet_name=4)
df=df[['KEGG.COMPOUND','CHEBI','LIPIDMAPS','PUBCHEM.COMPOUND','INCHIKEY','INCHIKEY.1']]

u = defaultdict(list)
for x in df:
    y =len(df[(df[x] == 1.0) | (df[x] == 2.0) |  (df[x] == 3.0)])/len(df)*100
    z =len(df[(df[x] == 1.0) | (df[x] == 2.0)])/len(df)*100
    w =len(df[df[x] == 1.0])/len(df)*100
    u[x].extend([y,z,w])

u ={k: v for k, v in sorted(u.items(), key=lambda item: item[1])}

l = []
l2 = [] 
for k, v in u.items():
    l.append(k)
    l2.append(v)
r = np.array(l2)

size=20
f, ax = plt.subplots(figsize=(20,9))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})
x=plt.barh(range(len(l)), r[:,0],color='#4A708B',tick_label=l)##1D2F6F#4A708B
y=plt.barh(range(len(l)), r[:,1],color='#8470FF',tick_label=l)
z=plt.barh(range(len(l)), r[:,2],color='#7A378B',tick_label=l)##8390FA
ax.tick_params(axis='both',labelsize=size)
# ax.legend([x,y,z], ['add (3)','replace (2)','no change (1)'],fontsize=size, title='Group',title_fontsize=size)
ax.legend([x,y,z], ['add','replace','no change'],fontsize=size, title='Group',title_fontsize=size)
plt.xlabel('xlabel', fontsize=size)
#plt.show()
plt.savefig("figures/metabolites_annot_group-THG-beta-1_vs_Human1.svg", format="svg")


# Compare Reactions: MB reactions, GPR, KEGG IDs (model1 vs model12) 

df=pd.read_excel(file,sheet_name=5)
df=df[['MASS_BALANCE GROUP [1,2]','KEGG.REACTION GROUP [1,2,3]','GENE ASSOCIATION GROUP [1,2,3]']]

u = defaultdict(list)
for x in df:
    y =len(df[(df[x] == 1.0) | (df[x] == 2.0) |  (df[x] == 3.0)])/len(df)*100
    z =len(df[(df[x] == 1.0) | (df[x] == 2.0)])/len(df)*100
    w =len(df[df[x] == 1.0])/len(df)*100
    u[x].extend([y,z,w])

u ={k: v for k, v in sorted(u.items(), key=lambda item: item[1])}

l = []
l2 = [] 
for k, v in u.items():
    l.append(k)
    l2.append(v)
r = np.array(l2)

size=20
f, ax = plt.subplots(figsize=(20,9))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})
x=plt.barh(range(len(l)), r[:,0],color='#4A708B',tick_label=l)##1D2F6F
y=plt.barh(range(len(l)), r[:,1],color='#8470FF',tick_label=l)
z=plt.barh(range(len(l)), r[:,2],color='#7A378B',tick_label=l)##8390FA
ax.tick_params(axis='both',labelsize=size)
#ax.legend([x,y,z], ['add (3)','replace (2)','no change (1)'],fontsize=size, title='Group',title_fontsize=size)
ax.legend([x,y,z], ['add','replace','no change'],fontsize=size, title='Group',title_fontsize=size)
plt.xlabel('xlabel', fontsize=size)
#plt.show()
plt.savefig("figures/reaction_annot_group-THG-beta-1_vs_Human1.svg", format="svg")


# Compare Genes: Ensembl IDs (model1 vs model12) 

df=pd.read_excel(file,sheet_name=6)
df=df[['ENSEMBLE GROUP [1,2,3]']]

u = defaultdict(list)
for x in df:
    y =len(df[(df[x] == 1.0) | (df[x] == 2.0) |  (df[x] == 3.0)])/len(df)*100
    z =len(df[(df[x] == 1.0) | (df[x] == 2.0)])/len(df)*100
    w =len(df[df[x] == 1.0])/len(df)*100
    u[x].extend([y,z,w])

u ={k: v for k, v in sorted(u.items(), key=lambda item: item[1])}

l = []
l2 = [] 
for k, v in u.items():
    l.append(k)
    l2.append(v)
r = np.array(l2)

size=20
f, ax = plt.subplots(figsize=(20,9))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})
x=plt.barh(range(len(l)), r[:,0],color='#4A708B',tick_label=l)##1D2F6F
y=plt.barh(range(len(l)), r[:,1],color='#8470FF',tick_label=l)
z=plt.barh(range(len(l)), r[:,2],color='#7A378B',tick_label=l)##8390FA
ax.tick_params(axis='both',labelsize=size)
#ax.legend([x,y,z], ['add (3)','replace (2)','no change (1)'],fontsize=size, title='Group',title_fontsize=size)
ax.legend([x,y,z], ['add','replace','no change'],fontsize=size, title='Group',title_fontsize=size)
plt.xlabel('xlabel', fontsize=size)
#plt.show()
plt.savefig("figures/genes_annot_group.svg", format="svg")


# Compare non annotated metabolites (model1 vs model12 vs model15 vs modelTHG)

f, ax = plt.subplots(figsize=(10,2))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})

l = ['%']
l2=[len([x for x in model1.metabolites if sorted([*x.annotation]) == ['bigg.metabolite', 'sbo', 'vmhmetabolite']])/len(model1.metabolites)*100,len([x for x in model12.metabolites if sorted([*x.annotation]) == ['bigg.metabolite', 'sbo', 'vmhmetabolite']])/len(model12.metabolites)*100,len([x for x in model15.metabolites if sorted([*x.annotation]) == ['bigg.metabolite', 'sbo', 'vmhmetabolite']])/len(model15.metabolites)*100]
r = np.array(l2)

size=15
x = plt.barh(range(len(l)), r[0],color='darkslateblue',tick_label=l)
y = plt.barh(range(len(l)), r[1],color='blueviolet',tick_label=l)
z = plt.barh(range(len(l)), r[2],color='darkviolet',tick_label=l)
ax.legend([z,y,x], ['THG-beta-2','THG-beta-1','Reference GEM'],fontsize=size, title='Model',title_fontsize=size)
plt.title('Metabolites with no metabolic database annotation', fontsize=size)
ax.tick_params(axis='both',labelsize=size)
plt.xticks([0,5,15,20,25,30])
plt.xlabel('', fontsize=size)
#plt.show()
plt.savefig("figures/non_annotated_metabolites-THG-beta-1_vs_THG-beta-1_vs_Human1.svg", format="svg")


# MEMOTE (model1 vs model12 vs model15 vs modelTHG)

memote=np.array([[ 81, 82,81,82], # total
       [ 73, 78,78,79],  # met
       [ 72, 72,72,72],  # reac 
       [46, 46,46,46], # gene 
       [87,90,88,91]]) # consistency

f, ax = plt.subplots(figsize=(20,9))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})

l= ['Total','Met annotation','Reac annotation','Gene annotation','Consistency']

size=17
h=plt.barh(range(len(l)), memote[0:,3],color='plum',tick_label=l)##8390FA
x=plt.barh(range(len(l)), memote[0:,2],color='darkviolet',tick_label=l)
y=plt.barh(range(len(l)), memote[0:,1],color='blueviolet',tick_label=l)
z=plt.barh(range(len(l)), memote[0:,0],color='darkslateblue',tick_label=l)##8390FA

ax.tick_params(axis='both',labelsize=size)
ax.legend([h,x,y,z], ['THG','THG-beta-2','THG-beta-1','Reference GEM'],fontsize=size, title='Model',title_fontsize=size)
plt.xlabel('%', fontsize=size)
plt.xticks([0,20,40,60,80,100])
plt.title('MEMOTE (categories)',size=17)
#plt.show()
plt.savefig("figures/MEMOTE-THG_vs_THG-beta-1_vs_THG-beta-1_vs_Human1.svg", format="svg")


# Compare components in model (model1 vs model12 vs model15 vs modelTHG)

x=np.array([len(model1.genes),len([x for x in model1.reactions if x.gene_reaction_rule]),len(model1.metabolites),len(model1.reactions)])
y=np.array([len(model12.genes),len([x for x in model12.reactions if x.gene_reaction_rule]),len(model12.metabolites),len(model12.reactions)])
z=np.array([len(model15.genes),len([x for x in model15.reactions if x.gene_reaction_rule]),len(model15.metabolites),len(model15.reactions)])
h=np.array([len(modelTHG.genes),len([x for x in modelTHG.reactions if x.gene_reaction_rule]),len(modelTHG.metabolites),len(modelTHG.reactions)])

l = ['Genes','Gene reaction associations','Metabolites','Reactions']
r = np.array([x,y,z,h]).T

size=20
f, ax = plt.subplots(figsize=(20,9))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})

x=plt.barh(range(len(l)), r[0:,3],color='plum',tick_label=l)
y=plt.barh(range(len(l)), r[0:,2],color='darkviolet',tick_label=l)
z=plt.barh(range(len(l)), r[0:,1],color='blueviolet',tick_label=l)##8390FA
h=plt.barh(range(len(l)), r[0:,0],color='darkslateblue',tick_label=l)##8390FA

ax.tick_params(axis='both',labelsize=size)
ax.legend([x,y,z,h], ['THG','THG-beta-2','THG-beta-1','Reference GEM'],fontsize=size, title='Model',title_fontsize=size)
plt.xlabel('', fontsize=size)
plt.title('About (model components)',fontsize=size)
#plt.show()
plt.savefig("figures/odel_components-THG_vs_THG-beta-1_vs_THG-beta-1_vs_Human1.svg", format="svg")


# Compare compartment in model (model12 vs model15)

y= defaultdict(list)

l = [x.id for x in model1.reactions]
l2 = list()

for x in model15.reactions:
    if not x.id in l:
        l2.append(x.id)
        y[str(x.annotation)].append(x.id)

z = [len(y[x]) for x in y]
z=Counter(z)

f, ax = plt.subplots(figsize=(17,5))
sns.set_theme(style="white",rc={'axes.spines.top': False,'axes.spines.right': False,'font.sans-serif':'DejaVu Sans'})

sns.barplot(x=list(z.keys()), y=list(z.values()), color="darkviolet")

size=18
ax.bar_label(ax.containers[0], label_type='edge',fontsize=size)
    
    
plt.ylabel('Unique reactions', fontsize=size)
plt.xlabel('THG-beta-1 vs THG-beta-1 Replicates', fontsize=size)
plt.yticks(fontsize=size)
plt.xticks(fontsize=size)
plt.yticks([0,300,600,900,1200,1500])
#plt.show()
plt.savefig("figures/isoenzyme_in_new_compartment_reactions-THG-beta-1_vs_THG-beta-2.svg", format="svg")


# Results of text similarity algorithm to identify metabolites

counts = [1104,313,0]
counts = pd.Series(counts,
                   index=[str(counts[0])+" RP", str(counts[1])+" WP", str(counts[2])+" NP"])
pie(counts,'Text similarity algorithm for metabolite identification (0.82)')


# Results of text JI-based algorithm to identify reactions

counts = [84.10,15.90]
counts = pd.Series(counts,
                   index=[str(counts[0])+" RP", str(counts[1])+" WP"])
pie2(counts,'JI-based algorithm for reaction identification')


# Results of gene-to-reaction algorithm to build SGPR and GPR

counts = [97.77,2.23,0]
counts = pd.Series(counts,
                   index=[str(counts[0])+" RP", str(counts[1])+" WP", str(counts[2])+" NP"])
pie(counts,'Algorithm-based automatic construction of gene-reaction associations')


# Results of text mass balance algorithm to mass balance metabolic reactions

counts = [82.1,13.85,0.85,0.24,2.96]
counts = pd.Series(counts,
                   index=[str(counts[0])+" Mass balanced", str(counts[1])+" a", str(counts[2])+" b", str(counts[3])+" c", str(counts[4])+" d"])
pie5(counts,'Algorithm-based automatic reaction mass balance')
