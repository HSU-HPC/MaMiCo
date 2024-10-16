import os

"""Utilities for interactive terminal applications."""


def select(
    options: list,
    pre_selected: int = None,
    title: str = None,
    clear_before: bool = True,
    clear_after: bool = True,
) -> int:
    """Displays a single item selection menu in the terminal.
    This function loops while no valid input has been received.
    Returns the index of the selected option or -1.

    Keyword arguemnts:
    options -- A list of options to select from (should be short and about the same length)
    pre_selected -- The index of the option to pre-select (confirm with enter) or None
    title -- The title displayed above the menu
    clear_before -- If true, the terminal will be cleared before displaying the menu
    clear_after -- If true, the terminal will be cleared after displaying the menu
    """
    prompt = f"\nSelect ({1}-{len(options)}): "
    if pre_selected is not None and pre_selected < 0:
        pre_selected += len(options)
    while True:
        if clear_before:
            clear()
            if title is not None:
                print(title)
            for i, option in enumerate(options):
                print(f"{ '>' if i == pre_selected else ' ' }[{i+1}]\t{option}")
        try:
            input_str = input(prompt).strip()
            if len(input_str) == 0 and pre_selected is not None:
                selected = pre_selected
            else:
                selected = int(input_str) - 1
            if selected >= 0 and selected < len(options):
                if clear_after:
                    clear()
                return selected
        except KeyboardInterrupt:
            if clear_after:
                clear()
            return -1
        except ValueError:
            pass


def clear() -> None:
    """Clears all output from the terminal."""
    os.system("cls" if os.name == "nt" else "clear")


def fmt_bold(text: str) -> None:
    """Returns text formatted to be printed to the terminal in bold font.
    On windows this has no effect.

    Keyword arguments:
    text -- The text to be formatted
    """
    if os.name == "nt":
        return text
    return f"\033[1m{text}\033[0m"


def fmt_red(text: str) -> None:
    """Returns text formatted to be printed to the terminal in red font.
    On windows this has no effect.

    Keyword arguments:
    text -- The text to be formatted
    """
    if os.name == "nt":
        return text
    return f"\033[0;31m{text}\033[0m"
