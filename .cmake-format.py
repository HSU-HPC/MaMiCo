# -----------------------------
# Options effecting formatting.
# -----------------------------
with section("format"):

  # How wide to allow formatted cmake files
  line_width = 120

  # How many spaces to tab for indent
  tab_size = 2

  # If true, separate flow control names from their parentheses with a space
  separate_ctrl_name_with_space = False

  # If true, separate function names from parentheses with a space
  separate_fn_name_with_space = False

  # If a statement is wrapped to more than one line, than dangle the closing
  # parenthesis on its own line.
  dangle_parens = False

  # If an argument group contains more than this many sub-groups (parg or kwarg
  # groups) then force it to a vertical layout.
  max_subgroups_hwrap = 6

  # If a positional argument group contains more than this many arguments, then
  # force it to a vertical layout.
  max_pargs_hwrap = 8

  # A list of command names which should always be wrapped
  always_wrap = ["add_executable", "add_custom_target", "FetchContent_Declare", "add_custom_command"]

