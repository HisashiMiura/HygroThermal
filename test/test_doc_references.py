"""コードのコメント・docstring の「式 (k)」と、計算方法の文書の評価方法の式の対応を確かめる。

式とコードの対応は、関数の docstring の Notes と、計算する行のコメントの「式 (k)」にだけ書き、対応の表は作らない
（設計の文書 2.7 節）。このテストは、次の 2 つを確かめる。

1. コードの「式 (k)」が、その節に実際にある式の番号か（番号が式の数を超えていないか）
2. 評価方法の各節のすべての式が、コードのどこかで参照されているか（_UNREFERENCED に挙げた式を除く）

「式 (k)」がどの節の式かは、次の順に決める。
- 同じ行で、「式 (k)」のあとに（「節の名前」）がある（例：式 (7)（「熱・水分の収支」））
- 同じ行で、「式 (k)」の前にいちばん近い「節の名前」がある（例：「湿り空気と水の性質」の式 (5)）
- モジュールの docstring に書いた、評価方法の節のファイル（docs/source/contents/03_xx_*.md）
「RR-110 の式 (5)」のように、節の名前でないものの「の式」は、文書の外の式として数えない。
番号の範囲の中でのずれ（式 (3) と式 (4) の取り違えなど）は、このテストでは見つけられない。
"""

import ast
import glob
import os
import re

import pytest

_ROOT = os.path.join(os.path.dirname(__file__), '..')
_CONTENTS = os.path.join(_ROOT, 'docs', 'source', 'contents')
_MODULES = os.path.join(_ROOT, 'modules')

# 評価方法の節の名前 → ファイル名（1 行目の見出しから作る）
_PAGES = {}
for _path in sorted(glob.glob(os.path.join(_CONTENTS, '03_*.md'))):
    with open(_path, encoding='utf-8') as _f:
        _PAGES[_f.readline().lstrip('#').strip()] = os.path.basename(_path)

# 式の並び：「式 (a)」に「・(b)」「〜(b)」が続く
_REF = re.compile(r'式 \((\d+)\)((?:[・〜]\(\d+\))*)')

# 節の名前（「」で囲んだもの）
_TITLE = re.compile(r'「([^」]+)」')

# 参照されていなくてよい式, {ファイル名: {式の番号}}（説明のための式など。足すときは理由を書く）
_UNREFERENCED: dict[str, set[int]] = {}


def _equation_count(file_name: str) -> int:
    """評価方法の節のファイルにある、番号の付いた式の数（「$$ (ラベル)」の数）。"""
    with open(os.path.join(_CONTENTS, file_name), encoding='utf-8') as f:
        return len(re.findall(r'^\$\$ \([^)]+\)', f.read(), flags=re.MULTILINE))


def _numbers(first: str, rest: str) -> list[int]:
    """「式 (a)・(b)〜(c)」の番号の並び（〜 はあいだの番号をすべて含める）。"""
    numbers = [int(first)]
    for sep, n in re.findall(r'([・〜])\((\d+)\)', rest):
        n = int(n)
        if sep == '〜':
            numbers.extend(range(numbers[-1] + 1, n + 1))
        else:
            numbers.append(n)
    return numbers


def _default_page(source: str) -> str | None:
    """モジュールの docstring に書いた、評価方法の節のファイル名（なければ None）。"""
    doc = ast.get_docstring(ast.parse(source)) or ''
    m = re.search(r'docs/source/contents/(03_\w+\.md)', doc)
    return m.group(1) if m else None


def _page_of(line: str, start: int, end: int, default: str | None) -> str | None:
    """行の start〜end にある「式 (k)」が、どの節の式か（文書の外の式なら ''、決まらなければ None）。"""

    # 「RR-110 の式」のように、節の名前でないものの「の式」は、文書の外の式
    before = line[:start]
    if before.endswith('の') and not before.endswith('」の'):
        return ''

    # 後ろの（「節の名前」）：次の「式 (」より前にあるもの
    after = line[end:]
    next_ref = after.find('式 (')
    after = after if next_ref < 0 else after[:next_ref]
    for m in _TITLE.finditer(after):
        if m.group(1) in _PAGES and after[:m.start()].rstrip().endswith('（'):
            return _PAGES[m.group(1)]

    # 前にいちばん近い「節の名前」
    titles = [m.group(1) for m in _TITLE.finditer(before) if m.group(1) in _PAGES]
    if titles:
        return _PAGES[titles[-1]]

    return default


def _references() -> list[tuple[str, int, str, int]]:
    """コードの「式 (k)」の一覧, [(モジュールのファイル名, 行, 節のファイル名, 式の番号)]。"""
    refs = []
    for path in sorted(glob.glob(os.path.join(_MODULES, '*.py'))):
        with open(path, encoding='utf-8') as f:
            source = f.read()
        default = _default_page(source)
        for i, line in enumerate(source.splitlines(), start=1):
            for m in _REF.finditer(line):
                page = _page_of(line, m.start(), m.end(), default)
                if page == '':
                    continue
                for n in _numbers(m.group(1), m.group(2)):
                    refs.append((os.path.basename(path), i, page, n))
    return refs


class TestParse:

    def test_numbers(self):
        # 「・」は並べ、「〜」はあいだの番号をすべて含める
        assert _numbers('1', '') == [1]
        assert _numbers('3', '・(4)') == [3, 4]
        assert _numbers('6', '〜(9)') == [6, 7, 8, 9]

    def test_page_of(self):
        # 後ろの（「節の名前」）、前の「節の名前」、モジュールの節、文書の外の式の順
        def page(line, default='03_05_eval_flows.md'):
            m = _REF.search(line)
            return _page_of(line, m.start(), m.end(), default)
        assert page('# 式 (7)（「熱・水分の収支」）') == '03_11_eval_balance.md'
        assert page('# 式 (9) の左辺 - 右辺（「熱・水分の収支」）') == '03_11_eval_balance.md'
        assert page('# 「湿り空気と水の性質」の式 (6)') == '03_02_eval_air_water.md'
        assert page('# 式 (9)') == '03_05_eval_flows.md'
        assert page('出典（RR-110 の式 (5)、ASHRAE）') == ''
        assert page('# 式 (9)', default=None) is None


def test_pages_found():
    # 評価方法の節が読み込めている（見出しの名前が変わると、ここで気づく）
    assert len(_PAGES) >= 13


def test_reference_page_is_known():
    # どの「式 (k)」も、どの節の式かが決まる
    unknown = [(module, line) for module, line, page, _ in _references() if page is None]
    assert unknown == [], f'どの節の式かが決まらない「式 (k)」：{unknown}'


def test_reference_number_in_range():
    # コードの「式 (k)」の番号が、その節の式の数を超えていない
    counts = {page: _equation_count(page) for page in _PAGES.values()}
    bad = [(module, line, page, n) for module, line, page, n in _references()
           if page is not None and not 1 <= n <= counts[page]]
    assert bad == [], f'節にない式の番号：{bad}'


@pytest.mark.parametrize('page', sorted(set(_PAGES.values()) - {'03_00_eval.md', '03_01_eval_flow.md'}))
def test_every_equation_referenced(page):
    # 評価方法の節のすべての式が、コードのどこかで参照されている
    referenced = {n for _, _, p, n in _references() if p == page}
    missing = set(range(1, _equation_count(page) + 1)) - referenced - _UNREFERENCED.get(page, set())
    assert missing == set(), f'{page} の、コードで参照されていない式：{sorted(missing)}'
