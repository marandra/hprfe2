#/usr/bin/env bash
_hprfe2_completions()
{
  COMPREPLY=($(compgen -W "init deploy generate pack resources" "${COMP_WORDS[1]}"))
}

complete -F _hprfe2_completions hprfe2
