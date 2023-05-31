# Set the shell to bash with printing of all commands (`-x`) and unofficial
# strict mode (`-euo pipefail`).
SHELL := bash -x -euo pipefail

# Default target, run on CI.  Runs all checks.
.PHONY: ci
ci: \
	lint-bash \
	check-format-bash \
	lint-awk \
	check-isort \
	check-black \
	check-snakefmt \
	flake8 \

# Alias lint to ci.
.PHONY: lint
lint: ci

# Run bash linting using spellcheck.
.PHONY: lint-bash
lint-bash:
	shellcheck -s bash scripts/*.sh

# Run bash formatting checks using beautysh.
.PHONY: check-format-bash
check-format-bash:
	for file in scripts/*.sh; do \
		beautysh --check --indent-size 4 $$file --force-function-style paronly; \
	done

# Run awk linting using gawk.
.PHONY: lint-awk
lint-awk:
	for file in scripts/*.awk; do \
		gawk -f $$file -Lfatal /dev/null >/dev/null; \
	done

# Run Python import sort checking with isort.
.PHONY: check-isort
check-isort:
	isort --profile=black --check-only tools scripts varfish_db_downloader

# Run Python format checking with black.
.PHONY: check-black
check-black:
	black --check --line-length 100 tools scripts varfish_db_downloader

# Run Snakemake format checking with snakefmt.
.PHONY: check-snakefmt
snakefmt:
	snakefmt --check --diff --line-length 100 Snakefile
	snakefmt --check --diff --line-length 100 snakefiles/*.smk

# Run Python linting with flake8.
.PHONY: flake8
flake8:
	flake8 --max-line-length 100 tools scripts varfish_db_downloader

# Run all automatic code formatting.
.PHONY: format
format:	\
	format-bash \
	isort \
	black \
	run-snakefmt

# Run bash formatting using beautysh.
.PHONY: format-bash
format-bash:
	for file in scripts/*.sh; do \
		beautysh --indent-size 4 $$file --force-function-style paronly; \
	done

# Run Python import sorting with isort.
.PHONY: isort
isort:
	isort --profile=black tools scripts varfish_db_downloader

# Run Python formatting with black.
.PHONY: black
black:
	black --line-length 100 tools scripts varfish_db_downloader

# Run Snakemake formatting with snakefmt.
.PHONY: run-snakefmt
run-snakefmt:
	snakefmt --line-length 100 Snakefile
	snakefmt --line-length 100 snakefiles/*.smk
