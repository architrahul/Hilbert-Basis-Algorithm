import requests
from bs4 import BeautifulSoup


def fetch_covering_online(v: int, k: int, t: int):
    url = f"https://ljcr.dmgordon.org/show_cover.php?v={v}&k={k}&t={t}"

    r = requests.get(
        url,
        headers={"User-Agent": "Mozilla/5.0"},
        timeout=30
    )
    r.raise_for_status()

    soup = BeautifulSoup(r.text, "lxml")
    pre = soup.find("pre")

    if pre is None:
        raise RuntimeError(f"No covering stored for C({v},{k},{t})")

    blocks = []
    for line in pre.text.strip().splitlines():
        row = [int(x) for x in line.split()]
        blocks.append(row)

    return blocks
