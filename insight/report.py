"""Page HTML « récit », autonome (images intégrées, aucun appel réseau) : les
figures dans l'ordre de lecture, chacune avec des constats chiffrés générés à
partir des résultats (pas d'interprétation inventée)."""
from __future__ import annotations

import base64
import datetime
import html
import os

CSS = """
:root{--bg:#fcfcfb;--card:#ffffff;--ink:#0b0b0b;--ink2:#52514e;--line:#e6e5e1;--acc:#2a78d6}
@media (prefers-color-scheme: dark){:root{--bg:#1a1a19;--card:#232321;--ink:#f0efec;
--ink2:#c3c2b7;--line:#383835;--acc:#3987e5} img{background:#fff;border-radius:6px}}
*{box-sizing:border-box}body{margin:0;background:var(--bg);color:var(--ink);
font:15px/1.55 system-ui,-apple-system,"Segoe UI",Roboto,sans-serif}
header{max-width:1100px;margin:0 auto;padding:40px 20px 8px}
h1{font-size:28px;margin:0 0 4px}header p{color:var(--ink2);margin:4px 0}
main{max-width:1100px;margin:0 auto;padding:8px 20px 60px}
section{background:var(--card);border:1px solid var(--line);border-radius:14px;
padding:22px 24px;margin:22px 0}h2{font-size:19px;margin:0 0 8px}
section img{width:100%;height:auto;margin-top:12px}ul{margin:8px 0 0 18px;padding:0}
li{margin:3px 0}.k{display:flex;flex-wrap:wrap;gap:12px;margin-top:14px}
.k div{border:1px solid var(--line);border-radius:10px;padding:10px 14px;min-width:150px}
.k b{display:block;font-size:22px}.k span{color:var(--ink2);font-size:13px}
nav{max-width:1100px;margin:0 auto;padding:0 20px;display:flex;flex-wrap:wrap;gap:8px}
nav a{color:var(--acc);text-decoration:none;font-size:14px;border:1px solid var(--line);
border-radius:20px;padding:4px 12px}footer{color:var(--ink2);font-size:13px;text-align:center}
em{font-style:italic}
"""


def _img(path):
    with open(path, "rb") as fh:
        return "data:image/png;base64," + base64.b64encode(fh.read()).decode()


def write_report(out_dir, title, kpis, sections, meta_lines):
    """sections : liste de dict(id, title, png, bullets)."""
    nav = "".join(f'<a href="#{s["id"]}">{html.escape(s["title"])}</a>' for s in sections)
    k = "".join(f"<div><b>{html.escape(str(v))}</b><span>{html.escape(n)}</span></div>"
                for n, v in kpis)
    body = []
    for s in sections:
        bl = "".join(f"<li>{b}</li>" for b in s.get("bullets", []))
        img = f'<img alt="{html.escape(s["title"])}" src="{_img(s["png"])}">' if s.get("png") else ""
        body.append(f'<section id="{s["id"]}"><h2>{html.escape(s["title"])}</h2>'
                    f'<ul>{bl}</ul>{img}</section>')
    meta = "".join(f"<p>{html.escape(m)}</p>" for m in meta_lines)
    page = (f"<!doctype html><html lang='en'><head><meta charset='utf-8'>"
            f"<meta name='viewport' content='width=device-width,initial-scale=1'>"
            f"<title>{html.escape(title)}</title><style>{CSS}</style></head><body>"
            f"<header><h1>{html.escape(title)}</h1>{meta}<div class='k'>{k}</div></header>"
            f"<nav>{nav}</nav><main>{''.join(body)}</main>"
            f"<footer>InSight · generated {datetime.date.today().isoformat()} · facts are "
            f"computed from the Deimos results; interpretation is yours</footer></body></html>")
    p = os.path.join(out_dir, "insight_report.html")
    with open(p, "w", encoding="utf-8") as fh:
        fh.write(page)
    return p
