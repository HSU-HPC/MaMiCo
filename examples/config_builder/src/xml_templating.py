from xml.dom import minidom


def remove_empty_files(s: str) -> str:
    return "\n".join([s for s in s.splitlines() if len(s.strip()) > 0])


def break_long_xml_lines(s: str) -> str:
    lines = s.splitlines()
    new_lines = []
    for line in lines:
        parts = line.split('" ')
        if len(parts) <= 2:
            new_lines.append(line)  # Avoid splitting trivial tags
            continue
        parts = [s + ('"' if i > 0 else "") for i, s in enumerate(parts[::-1])][::-1]
        new_lines.append(parts[0])
        # Indent like tag
        indentation = parts[0].split("<")[0]
        # Indent by tag name length
        tag_name = parts[0][len(indentation) :].split(" ")[0]
        prefix = indentation + " " * (len(tag_name) + 1)
        for part in parts[1:]:
            if part.endswith("/>"):
                new_lines.append(prefix + part[:-2])
                new_lines.append(indentation + "/>")
            else:
                new_lines.append(prefix + part)
    return "\n".join(new_lines)


def format_xml(s: str) -> str:
    s = minidom.parseString(s).toprettyxml(indent=" " * 4)
    s = remove_empty_files(s)
    s = break_long_xml_lines(s)
    return s


class PartialXml:
    def __init__(self, xml_template_str) -> None:
        self._string = xml_template_str

    def substitute(self, a: str, b: object) -> None:
        self._string = self._string.replace(a, str(b))

    def get(self) -> str:
        if "{" in self._string:
            raise RuntimeError("Not all values have been substituted!")
        return format_xml(self._string)
