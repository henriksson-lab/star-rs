#!/usr/bin/env python3
import json
import keyword
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
REPORT = Path("/tmp/star_cpp.json")
SRC = ROOT / "STAR" / "source"

RUST_KEYWORDS = {
    "as", "break", "const", "continue", "crate", "else", "enum", "extern",
    "false", "fn", "for", "if", "impl", "in", "let", "loop", "match", "mod",
    "move", "mut", "pub", "ref", "return", "self", "Self", "static",
    "struct", "super", "trait", "true", "type", "unsafe", "use", "where",
    "while", "async", "await", "dyn", "abstract", "become", "box", "do",
    "final", "macro", "override", "priv", "typeof", "unsized", "virtual",
    "yield", "try",
}


def snake(value: str) -> str:
    value = value.replace("::", "_")
    value = re.sub(r"[^A-Za-z0-9_]+", "_", value)
    value = re.sub(r"_+", "_", value).strip("_")
    if not value:
        value = "unnamed"
    if value[0].isdigit():
        value = "n_" + value
    value = value.lower()
    if value in RUST_KEYWORDS or keyword.iskeyword(value):
        value += "_"
    return value


def rust_type(text: str, category: str) -> str:
    t = text.replace(" ", "")
    if category == "bool" or t == "bool":
        return "bool"
    if category == "string" or "string" in t:
        return "String"
    if category == "float" or any(x in t for x in ("float", "double")):
        return "f64"
    if category == "pointer" or "*" in t:
        return "Option<usize>"
    if "uint64" in t or "uint64_t" in t or t in {"uint", "uint32", "uint32_t"}:
        return "u64" if "64" in t else "u32"
    if "int64" in t or "int64_t" in t:
        return "i64"
    if "int" in t or t in {"long", "short"}:
        return "i64"
    if category in {"array", "collection"} or any(x in t for x in ("vector<", "array<", "map<", "set<")):
        return "Vec<usize>"
    return "()"


def rust_type_name(name: str) -> str:
    rust_name = re.sub(r"[^A-Za-z0-9_]", "_", name)
    if not rust_name or rust_name[0].isdigit():
        rust_name = "Type_" + rust_name
    return rust_name


def discover_named_types(report):
    types = {}
    for item in report.get("structs", []):
        types[item["name"]] = item

    pattern = re.compile(r"^\s*(class|struct)\s+([A-Za-z_][A-Za-z0-9_]*)\b")
    for path in SRC.rglob("*"):
        if "htslib" in path.parts or not path.suffix in {".h", ".cpp", ".hpp"}:
            continue
        for lineno, line in enumerate(path.read_text(errors="ignore").splitlines(), 1):
            m = pattern.match(line)
            if not m:
                continue
            name = m.group(2)
            if name in {"stat", "statvfs", "shmid_ds"}:
                continue
            types.setdefault(name, {
                "name": name,
                "kind": m.group(1),
                "location": {"file": str(path.relative_to(ROOT)), "line_start": lineno},
                "fields": [],
            })
    return dict(sorted(types.items()))


def write_structs(types):
    out = ROOT / "src" / "generated" / "structs.rs"
    lines = [
        "#![allow(non_camel_case_types)]",
        "#![allow(non_snake_case)]",
        "",
    ]
    used_type_names = set()
    emitted = {}
    for name, item in types.items():
        rust_name = rust_type_name(name)
        base = rust_name
        i = 2
        while rust_name in used_type_names:
            rust_name = f"{base}_{i}"
            i += 1
        used_type_names.add(rust_name)
        emitted[name] = rust_name
        loc = item.get("location", {})
        lines.append(f"#[doc = \"Original {item.get('kind', 'type')} `{name}` at {loc.get('file', '')}:{loc.get('line_start', '')}.\"]")
        lines.append("#[derive(Clone, Debug, Default, PartialEq)]")
        lines.append(f"pub struct {rust_name} {{")
        used_fields = set()
        for field in item.get("fields", []):
            fname = snake(field.get("name", "field"))
            base_field = fname
            i = 2
            while fname in used_fields:
                fname = f"{base_field}_{i}"
                i += 1
            used_fields.add(fname)
            ty = rust_type(field.get("ty", {}).get("text", ""), field.get("category", "other"))
            lines.append(f"    pub {fname}: {ty},")
        lines.append("}")
        lines.append("")
    out.write_text("\n".join(lines))
    return emitted


def write_functions(functions, struct_names):
    out = ROOT / "src" / "generated" / "functions.rs"
    lines = [
        "#![allow(non_snake_case)]",
        "",
    ]
    mappings = [
        "# Generated from STAR/source with code-complexity-comparator.",
        "# Rust names include source file and line to disambiguate C++ overloads.",
        "",
    ]
    names = set()
    for f in functions:
        loc = f["location"]
        stem = snake(Path(loc["file"]).stem)
        base = f"{stem}_l{loc['line_start']}_{snake(f['name'])}"
        name = base
        i = 2
        while name in names:
            name = f"{base}_{i}"
            i += 1
        names.add(name)
        inputs = f.get("signature", {}).get("inputs", [])
        arg_note = ", ".join(
            f"{a.get('name', '')}: {a.get('ty', {}).get('text', '')}".strip()
            for a in inputs
        )
        lines.append(f"#[doc = \"Original `{f['name']}` at {loc['file']}:{loc['line_start']}. Args: {arg_note}\"]")
        lines.append(f"pub fn {name}() {{")
        lines.append(f"    todo!(\"translate {f['name']} from {loc['file']}:{loc['line_start']}\");")
        lines.append("}")
        lines.append("")
        mappings.extend([
            "[[entries]]",
            f"rust = \"{name}\"",
            "rust_path = \"src/generated/functions.rs\"",
            f"other = {json.dumps(f['name'])}",
            f"other_path = {json.dumps(loc['file'])}",
            f"other_line = {loc['line_start']}",
            "",
        ])
    report_structs = [name for name in struct_names if any(ch in name for ch in "<>:")]
    if report_structs:
        mappings.extend([
            "# Struct mappings for C++ names that are not legal Rust identifiers.",
            "",
        ])
    for other_name in report_structs:
        mappings.extend([
            "[[entries]]",
            f"rust = {json.dumps(struct_names[other_name])}",
            "rust_path = \"src/generated/structs.rs\"",
            f"other = {json.dumps(other_name)}",
            "",
        ])
    out.write_text("\n".join(lines))
    (ROOT / "ccc_mapping.toml").write_text("\n".join(mappings))


def main():
    report = json.loads(REPORT.read_text())
    (ROOT / "src" / "generated").mkdir(parents=True, exist_ok=True)
    (ROOT / "ccc").mkdir(exist_ok=True)
    (ROOT / "ccc" / "star_cpp.json").write_text(json.dumps(report, indent=2) + "\n")
    order = Path("/tmp/star_order.csv")
    if order.exists():
        (ROOT / "ccc" / "porting_order.csv").write_text(order.read_text())
    struct_names = write_structs(discover_named_types(report))
    write_functions(report["functions"], struct_names)


if __name__ == "__main__":
    main()
