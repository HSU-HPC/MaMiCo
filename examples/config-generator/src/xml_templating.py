from xml.dom import minidom

"""Utilities for working with XML (templates)."""


def remove_empty_files(s: str) -> str:
    """Removes all lines containing only whitespace from a string."""
    return "\n".join([s for s in s.splitlines() if len(s.strip()) > 0])


def break_long_xml_lines(xml_str: str) -> str:
    """Format a serialized XML document by putting every key in a separat line and neatly aligning them."""
    lines = xml_str.splitlines()
    new_lines = []
    for line in lines:
        if line.strip().startswith("<?"):
            new_lines.append(line)  # Avoid splitting meta-tags
            continue
        parts = line.split('" ')
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


def format_xml(xml_str: str) -> str:
    """Prettify a serialized XML document to be rather tall than wide."""
    xml_str = minidom.parseString(xml_str).toprettyxml(indent=" " * 4)
    xml_str = remove_empty_files(xml_str)
    xml_str = break_long_xml_lines(xml_str)
    return xml_str


class PartialXml:
    """A class to facilitate substitution of XML templates (using "{placeholder}")"""

    def __init__(self, xml_template_str) -> None:
        """Create a new instance from an XML template string (using "{placeholder}")"""
        self._string = xml_template_str

    def substitute(self, placeholder: str, value: object) -> None:
        """Create a new instance from an XML template string (using "{placeholder}")

        Keyword arguments:
        placeholder -- The concrete placeholder ("{placeholder}") to replace in the template
        value -- The concrete value to replace "{placeholder}" with in the template
        """
        self._string = self._string.replace("{" + placeholder + "}", str(value))

    def get(self) -> str:
        """Returns the current value of the template.
        A RuntimeError occurs if not all placeholders have been substituted."""
        if "{" in self._string:
            raise RuntimeError("Not all values have been substituted!")
        return format_xml(self._string)
