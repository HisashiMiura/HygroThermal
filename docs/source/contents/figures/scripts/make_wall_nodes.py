"""評価方法「壁の構成と質点」（03_04）の図 1（層と質点・つながり・コントロールボリューム、水膜の質点を含む）の SVG を作る。

使い方（リポジトリの一番上のフォルダで）：
    python docs/source/contents/figures/scripts/make_wall_nodes.py

figures/wall_nodes.svg を書き出す。図の大きさ（層・境目の帯・通気層の幅など）は、下の定数で変える。
つながりの抵抗は zigzag.py のギザギザで描く（幅一定、線の中央。線が短いときだけ縮める）。
"""
import os

from zigzag import zigzag

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'wall_nodes.svg')

ORANGE, GRAY, INK, INK2, FILM, AIR = '#d9822b', '#9a9893', '#2b2a28', '#52514e', '#1b8a6b', '#2f6f9f'
Y = 140          # つながりの線の高さ
BAND = 90        # 広げて描く境目の帯の幅（水膜の質点は帯の中央）
MAT = 220        # 材料の層の幅（分割数 3：質点の間隔 110）
AIRW = 200       # 通気層の幅
SURF = 70        # 外気・室内とのつながりの長さ

x0 = 100                       # 外表面
m0_r = x0 + MAT               # 層 0 の右の表面
f1 = m0_r + BAND / 2          # 水膜 1
a_l = m0_r + BAND             # 通気層の左
air = a_l + AIRW / 2          # 空気の質点
a_r = a_l + AIRW              # 通気層の右
f2 = a_r + BAND / 2           # 水膜 2
m2_l = a_r + BAND             # 層 2 の左の表面
x_in = m2_l + MAT             # 内表面
W = x_in + SURF + 40          # 図の幅（右の余白を含む）

nodes = [(x0, 'm'), (x0 + MAT / 2, 'm'), (m0_r, 'm'), (f1, 'f'), (air, 'a'), (f2, 'f'), (m2_l, 'm'), (m2_l + MAT / 2, 'm'), (x_in, 'm')]
links = [(x0 - SURF, x0, ORANGE)] + \
        [(nodes[i][0], nodes[i + 1][0], GRAY if nodes[i][1] == 'm' and nodes[i + 1][1] == 'm' else ORANGE) for i in range(8)] + \
        [(x_in, x_in + SURF, ORANGE)]

o = []
w = o.append
w(f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {W:g} 350" width="{W:g}" height="350" '
  f'font-family="\'Hiragino Sans\', \'Noto Sans JP\', \'Yu Gothic\', sans-serif" font-size="13" fill="{INK}">')
w('  <!-- 例：層 0 = 材料（分割数 3）、層 1 = 通気層、層 2 = 材料（分割数 3）。層と層の境目には水膜の質点を置く。'
  '質点は N = 9 個、つながりは 10 本。境目は厚さがないが、水膜の質点が見えるように、境目を広げて（帯で）描く。'
  'つながりの抵抗は、幅一定のギザギザを線の中央に描く（線が短いときだけ縮める） -->')
# 背景
w(f'  <rect x="{x0}" y="70" width="{MAT}" height="200" fill="#f3ece0"/>')
w(f'  <rect x="{m0_r}" y="70" width="{BAND}" height="200" fill="#e3f1ec"/>')
w(f'  <rect x="{a_l}" y="70" width="{AIRW}" height="200" fill="#e6eff6"/>')
w(f'  <rect x="{a_r}" y="70" width="{BAND}" height="200" fill="#e3f1ec"/>')
w(f'  <rect x="{m2_l}" y="70" width="{MAT}" height="200" fill="#f3ece0"/>')
# 層の名前
for cx, t in ((x0 + MAT / 2, '0（材料）'), (air, '1（通気層）'), (m2_l + MAT / 2, '2（材料）')):
    w(f'  <text x="{cx:g}" y="58" text-anchor="middle">層 <tspan font-style="italic">l</tspan> = {t}</text>')
# 外気・室内
w(f'  <text x="{x0 - SURF + 10:g}" y="162" text-anchor="middle">外気</text>')
w(f'  <text x="{x_in + SURF - 10:g}" y="162" text-anchor="middle">室内</text>')
# 表面の線
w(f'  <g stroke="{INK}" stroke-width="2.5">')
for x in (x0, m0_r, a_l, a_r, m2_l, x_in):
    w(f'    <line x1="{x:g}" y1="70" x2="{x:g}" y2="270"/>')
w('  </g>')
w(f'  <g text-anchor="middle" font-size="12" fill="{INK2}">')
for x, t in ((x0, '外表面'), (f1, '境目'), (f2, '境目'), (x_in, '内表面')):
    w(f'    <text x="{x:g}" y="290">{t}</text>')
w('  </g>')
# Δx_0
dx_r = x0 + MAT / 2
w(f'  <g stroke="{INK2}" stroke-width="1" fill="none">')
w(f'    <line x1="{x0}" y1="96" x2="{dx_r:g}" y2="96"/>')
w(f'    <path d="M{x0 + 6},92 L{x0},96 L{x0 + 6},100 M{dx_r - 6:g},92 L{dx_r:g},96 L{dx_r - 6:g},100"/>')
w('  </g>')
w(f'  <text x="{(x0 + dx_r) / 2:g}" y="90" text-anchor="middle" font-size="12">Δ<tspan font-style="italic">x</tspan><tspan font-size="9" dy="3">0</tspan></text>')
# つながり
for a, b, col in links:
    w('  ' + zigzag(a, b, Y, col))
# つながりの番号
w('  <g text-anchor="middle" font-size="12">')
for k, (a, b, _) in enumerate(links):
    w(f'    <text x="{(a + b) / 2:g}" y="124"><tspan font-style="italic">k</tspan>={k}</text>')
w('  </g>')
# 質点
w(f'  <g stroke="#fcfcfb" stroke-width="2">')
for x, t in nodes:
    if t == 'f':
        w(f'    <path d="M{x:g},{Y - 8} L{x + 8:g},{Y} L{x:g},{Y + 8} L{x - 8:g},{Y} Z" fill="{FILM}"/>')
    else:
        w(f'    <circle cx="{x:g}" cy="{Y}" r="6" fill="{AIR if t == "a" else INK}"/>')
w('  </g>')
# 質点の番号
w('  <g fill="#fcfcfb">')
for i, (x, t) in enumerate(nodes):
    y = 172 if t == 'f' else 154
    w(f'    <rect x="{x - 17:g}" y="{y}" width="34" height="16" rx="3"/>')
w('  </g>')
w('  <g text-anchor="middle" font-size="12">')
for i, (x, t) in enumerate(nodes):
    y = 184 if t == 'f' else 166
    fill = f' fill="{FILM}"' if t == 'f' else ''
    w(f'    <text x="{x:g}" y="{y}"{fill}><tspan font-style="italic">i</tspan>={i}</text>')
w('  </g>')
# コントロールボリューム
cvs = [(x0, x0 + MAT / 4), (x0 + MAT / 4, x0 + 3 * MAT / 4), (x0 + 3 * MAT / 4, m0_r), None, (a_l, a_r), None,
       (m2_l, m2_l + MAT / 4), (m2_l + MAT / 4, m2_l + 3 * MAT / 4), (m2_l + 3 * MAT / 4, x_in)]
w(f'  <g fill="none" stroke="{INK2}" stroke-width="1" stroke-dasharray="4 3">')
for cv in cvs:
    if cv:
        w(f'    <rect x="{cv[0] + 2:g}" y="198" width="{cv[1] - cv[0] - 4:g}" height="26"/>')
w('  </g>')
w(f'  <g stroke="{FILM}" stroke-width="1.5">')
for x in (f1, f2):
    w(f'    <line x1="{x:g}" y1="198" x2="{x:g}" y2="224"/>')
w('  </g>')
w(f'  <g text-anchor="middle" font-size="11" fill="{INK2}">')
for i, cv in enumerate(cvs):
    if cv:
        w(f'    <text x="{(cv[0] + cv[1]) / 2:g}" y="216"><tspan font-style="italic">w</tspan><tspan font-size="8" dy="3">{i}</tspan></text>')
w('  </g>')
w(f'  <g text-anchor="middle" font-size="10" fill="{FILM}">')
for i, x in ((3, f1), (5, f2)):
    w(f'    <text x="{x:g}" y="238"><tspan font-style="italic">w</tspan><tspan font-size="8" dy="3">{i}</tspan><tspan dy="-3"> = 0</tspan></text>')
w('  </g>')
# 質点 1 の d_L,1 と d_R,1
n1 = x0 + MAT / 2; l1 = x0 + MAT / 4; r1 = x0 + 3 * MAT / 4
w(f'  <g stroke="{INK2}" stroke-width="1" fill="none">')
w(f'    <line x1="{l1:g}" y1="250" x2="{n1:g}" y2="250"/>')
w(f'    <line x1="{n1:g}" y1="250" x2="{r1:g}" y2="250"/>')
w(f'    <path d="M{l1 + 6:g},246 L{l1:g},250 L{l1 + 6:g},254 M{n1 - 6:g},246 L{n1:g},250 L{n1 - 6:g},254 '
  f'M{n1 + 6:g},246 L{n1:g},250 L{n1 + 6:g},254 M{r1 - 6:g},246 L{r1:g},250 L{r1 - 6:g},254"/>')
w(f'    <line x1="{n1:g}" y1="230" x2="{n1:g}" y2="254" stroke-dasharray="2 2"/>')
w('  </g>')
w('  <g text-anchor="middle" font-size="11">')
w(f'    <text x="{(l1 + n1) / 2:g}" y="266"><tspan font-style="italic">d</tspan><tspan font-size="8" dy="3">L,1</tspan></text>')
w(f'    <text x="{(n1 + r1) / 2:g}" y="266"><tspan font-style="italic">d</tspan><tspan font-size="8" dy="3">R,1</tspan></text>')
w('  </g>')
# 凡例
lx = 60
w('  <g font-size="12">')
w('    ' + zigzag(lx, lx + 60, 324, ORANGE)); w(f'    <text x="{lx + 68}" y="328">境目の抵抗が入るつながり</text>')
lx = 310
w('    ' + zigzag(lx, lx + 60, 324, GRAY)); w(f'    <text x="{lx + 68}" y="328">層の中のつながり</text>')
lx = 510
w(f'    <rect x="{lx}" y="316" width="30" height="16" fill="none" stroke="{INK2}" stroke-dasharray="4 3"/>')
w(f'    <text x="{lx + 38}" y="328">コントロールボリューム</text>')
lx = 760
w(f'    <path d="M{lx},316 L{lx + 8},324 L{lx},332 L{lx - 8},324 Z" fill="{FILM}"/>')
w(f'    <text x="{lx + 16}" y="328">水膜の質点（幅 0）</text>')
w('  </g>')
w('</svg>')
open(OUT, 'w').write('\n'.join(o) + '\n')
print(f'{os.path.normpath(OUT)} を書き出しました（幅 {W:g} px）。')
