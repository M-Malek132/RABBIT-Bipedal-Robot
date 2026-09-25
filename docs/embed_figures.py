#!/usr/bin/env python3
"""Re-embed report figures into the self-contained HTML reports.

    python3 docs/embed_figures.py docs/ch4_report.html
    python3 docs/embed_figures.py --check docs/ch4_report.html

The English reports carry their figures as base64 data URIs so that one .html
file is the whole report. That makes them easy to send and impossible to keep
current by hand: every time ch4_doc_figures redraws docs/figures/*.png, each
copy inside the HTML has to be replaced too, and nothing says when one was
missed.

Every <img> that names its source with data-fig="figures/<name>.png" (the path
is relative to the HTML file) gets its src rebuilt from that file. Images
without data-fig are left alone. --check writes nothing and exits 1 when any
embedded copy differs from its file, so it can be run after the figure scripts
to see whether the report is stale.

Standard library only.
"""

import argparse
import base64
import os
import re
import sys
import tempfile

MIME = {
    '.png': 'image/png',
    '.jpg': 'image/jpeg',
    '.jpeg': 'image/jpeg',
    '.gif': 'image/gif',
    '.svg': 'image/svg+xml',
}

IMG = re.compile(r'<img\b[^>]*>', re.IGNORECASE)
FIG = re.compile(r'\bdata-fig="([^"]+)"')
SRC = re.compile(r'\bsrc="[^"]*"')


def data_uri(path):
    ext = os.path.splitext(path)[1].lower()
    if ext not in MIME:
        raise ValueError(f'{path}: no MIME type for "{ext}"')
    with open(path, 'rb') as fh:
        payload = base64.b64encode(fh.read()).decode('ascii')
    return f'data:{MIME[ext]};base64,{payload}'


def process(html_path, check):
    base = os.path.dirname(os.path.abspath(html_path))
    with open(html_path, encoding='utf-8') as fh:
        text = fh.read()

    stats = {'updated': [], 'current': [], 'missing': []}

    def rebuild(match):
        tag = match.group(0)
        fig = FIG.search(tag)
        if not fig:
            return tag
        rel = fig.group(1)
        path = os.path.join(base, rel)
        if not os.path.isfile(path):
            stats['missing'].append(rel)
            return tag
        new_src = f'src="{data_uri(path)}"'
        old_src = SRC.search(tag)
        if old_src is None:
            new_tag = tag.replace('<img', f'<img {new_src}', 1)
        else:
            new_tag = tag[:old_src.start()] + new_src + tag[old_src.end():]
        stats['current' if new_tag == tag else 'updated'].append(rel)
        return new_tag

    new_text = IMG.sub(rebuild, text)

    name = os.path.basename(html_path)
    n = sum(len(v) for v in stats.values())
    print(f'{name}: {n} figure(s) named by data-fig')
    for rel in stats['updated']:
        print(f'  {"STALE  " if check else "updated"} {rel}')
    for rel in stats['current']:
        print(f'  current {rel}')
    for rel in stats['missing']:
        print(f'  MISSING {rel}')

    if not check and new_text != text:
        # Write beside the target and rename over it, so an interrupted run
        # never leaves a half-written report.
        fd, tmp = tempfile.mkstemp(dir=base, suffix='.html.tmp')
        with os.fdopen(fd, 'w', encoding='utf-8') as fh:
            fh.write(new_text)
        os.replace(tmp, html_path)

    return stats


def main():
    ap = argparse.ArgumentParser(description=__doc__.split('\n\n')[0])
    ap.add_argument('html', nargs='+', help='report(s) to update')
    ap.add_argument('--check', action='store_true',
                    help='report stale figures without writing; exit 1 if any')
    args = ap.parse_args()

    bad = False
    for html in args.html:
        stats = process(html, args.check)
        if stats['missing'] or (args.check and stats['updated']):
            bad = True
    sys.exit(1 if bad else 0)


if __name__ == '__main__':
    main()
