# Set the shell to bash with printing of all commands (`-x`) and unofficial
# strict mode (`-euo pipefail`).
SHELL := bash -x -euo pipefail

# Default target, run on CI.  Runs all checks.
.PHONY: ci
ci: \
	lint-bash \
	check-format-bash \
	lint-awk

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

# Run all automatic code formatting.
.PHONY: format
format:	\
	format-bash

# Run bash formatting using beautysh.
.PHONY: format-bash
format-bash:
	for file in scripts/*.sh; do \
		beautysh --indent-size 4 $$file --force-function-style paronly; \
	done
