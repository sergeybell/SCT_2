# Contributing to SCT_2

Workflow for a change to this repository. Agent rules are in [`AGENTS.md`](AGENTS.md).

## Workflow

One issue, one branch, one pull request.

1. Open an issue, or take an open one.
2. Branch from current `main`. The name is `<type>/<short-description>`, for example `feat/magnetic-ring-lattice`.
3. Open a pull request that closes that issue (`Closes #N`).
4. Short acceptance checks from the issue are green (script run, numbers table, figures).
5. Squash-merge into `main`.

`main` receives squash-merges only. Do not force-push, and do not rewrite the history of `main`.

## Conventional Commits

Commits follow [Conventional Commits 1.0.0](https://www.conventionalcommits.org/en/v1.0.0/).

```text
<type>[optional scope][optional !]: <description>
```

Types: `feat`, `fix`, `refactor`, `docs`, `test`, `chore`, `perf`, `ci`, `build`.

A breaking change is marked with `!` after the type (and the scope, if there is one), or with a `BREAKING CHANGE:` footer.

Squash-merge keeps one commit on `main`, and that commit is the pull request title. The title is a conventional commit. The body says why.

## Checks

There is no package-wide test suite for SCT_2. A study PR is accepted when:

- the step script runs with the short default settings;
- the numbers table and short-run figures are in the PR;
- any BELL findings are listed (draft text for a BELL issue, not a patch there).

Long tracking runs are behind an explicit flag. Progress must be printed unbuffered (`flush`).

## Changelog

This repository does not keep a changelog. Do not add `CHANGELOG.md`.
