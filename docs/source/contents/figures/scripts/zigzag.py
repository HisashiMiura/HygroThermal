"""抵抗の記号（ギザギザ）の SVG の path を作る（計算方法の文書の図で使う）。

ギザギザの部分の幅は一定（ZW）で線の中央に置き、線が短くて入らないときだけ、両側に LEAD を残して縮める。
"""
ZW = 40.0     # ギザギザの部分の幅, px
HALF = 6      # 山と谷の数（半波の数）
AMP = 5.0     # 振れ幅, px
LEAD = 6.0    # ギザギザの両側に、少なくとも残す直線, px

def zigzag(x1: float, x2: float, y: float, color: str, width: float = 2) -> str:
    """x1 から x2 までの水平な線（高さ y）を、中央にギザギザを入れた SVG の path にする。"""
    L = x2 - x1
    zw = min(ZW, L - 2 * LEAD)
    a = x1 + (L - zw) / 2
    step = zw / HALF
    pts = [f"M{x1:g},{y:g}", f"L{a:.1f},{y:g}"]
    for k in range(HALF):
        x = a + step * (k + 0.5)
        yy = y - AMP if k % 2 == 0 else y + AMP
        pts.append(f"L{x:.1f},{yy:g}")
    pts.append(f"L{a + zw:.1f},{y:g}")
    pts.append(f"L{x2:g},{y:g}")
    return (f'<path d="{" ".join(pts)}" fill="none" stroke="{color}" stroke-width="{width}" '
            f'stroke-linejoin="round" stroke-linecap="round"/>')
