"""入力の dict を確かめる共通の部品（壁の入力と、計算の入力で使う）。"""


def check_keys(d: dict, required: tuple[str, ...], where: str, optional: tuple[str, ...] = ()):
    """dict のキーを確かめる。必ず指定するキーがない、または知らないキー（打ち間違いなど）があればエラーにする。

    Args:
        d: 入力の dict
        required: 必ず指定するキー
        where: エラーのメッセージに出す、どの入力か
        optional: 省いてよいキー（必ず指定するキーと、省いてよいキーのほかは使えない）
    """

    # 必ず指定するキーのうち、dict にないものを集める。1 つでもあればエラー
    missing = [k for k in required if k not in d]
    if missing:
        raise ValueError(f'{where} に {missing} がありません。')

    # dict のキーのうち、使えるキーにないもの（打ち間違いなど）を集める。1 つでもあればエラー
    unknown = sorted(set(d) - set(required) - set(optional))
    if unknown:
        raise ValueError(f'{where} に知らないキー {unknown} があります（使えるキー: {list(required) + list(optional)}）。')


def check_number(value, name: str, where: str) -> float:
    """数（int か float。True・False は数として扱わない）であることを確かめて、float にして返す。

    Args:
        value: 入力の値
        name: エラーのメッセージに出す、キーの名前
        where: エラーのメッセージに出す、どの入力か
    """

    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ValueError(f'{where} の {name} は数にしてください（{value!r}）。')
    return float(value)
