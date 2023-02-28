# Memote tasks

This repository is a python module `metabolic_tasks` and extra tests for memote `tests_extra` that will use the python module to run tests perfoming metabolic tasks when called from memote.

## Memote setup

Usually, you would need to start from a memote repository, so the first step is to start memote. Here, while [this PR](https://github.com/opencobra/memote/pull/729) is not merged, I would recomend to install it from the fork:

```bash
pip install "git+https://github.com/carrascomj/memote.git@refactor-github-actions"
```

Then run

```bash
memote new
```

As usual.

## Running the tasks standlone

To run the whole test suite with metabolic tasks:

```bash
memote report snapshot --custom-tests tests_extra
```

To run only the metabolic tasks:

```bash
memote report snapshot --custom-tests tests_extra --exclusive test_essential_metabolic_tasks
```

## Integrate the tasks into memote CI/CD

The easiest is to copy and push both folders to `gh-pages` (deployment branch) once it has been pushed once (ran `memote online`).

Afterwards, you would need to modify [this line](https://github.com/carrascomj/cookiecutter-memote/blob/refactor-gh-actions/%7B%7Bcookiecutter.project_slug%7D%7D/scripts/ghactions_deploy.sh#L34) in your local repository at `scripts/ghactions_deploy.sh` and push to master:

```diff
- memote report history --filename="${output}"
+ memote report history --custom-tests tests_extra --filename="${output}"
```