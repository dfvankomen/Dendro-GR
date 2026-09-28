# Option-group helpers shared by the solvers.

# Fails if more than one option in the group is ON. Use for groups the source
# resolves with an #ifdef/#elif chain, where a second option is silently ignored.
function(dgr_require_at_most_one LABEL)
  set(_on "")
  foreach(_opt ${ARGN})
    if(${_opt})
      list(APPEND _on "${_opt}")
    endif()
  endforeach()
  list(LENGTH _on _count)
  if(_count GREATER 1)
    string(REPLACE ";" "\n  " _pretty "${_on}")
    message(FATAL_ERROR
            "${LABEL}: enable at most one of these, got ${_count}:\n  ${_pretty}")
  endif()
endfunction()
