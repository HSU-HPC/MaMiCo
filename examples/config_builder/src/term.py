import os


def select(
    options: list,
    pre_selected: int = None,
    title: str = None,
    clear_before: bool = True,
    clear_after: bool = True,
) -> int:
    prompt = f"\nSelect ({1}-{len(options)}): "
    if pre_selected < 0:
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
    os.system("cls" if os.name == "nt" else "clear")
