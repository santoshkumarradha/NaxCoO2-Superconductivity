#!/usr/bin/env python3
"""Minimal fal.ai gemini-3-pro-image-preview generator (stdlib only).

Usage:
  _falgen.py --image COND.png --prompt-file P.txt --out OUT.png \
             [--aspect 21:9] [--resolution 4K] [--seed N]
"""
import argparse, base64, json, mimetypes, os, sys, urllib.request, urllib.error

ENDPOINT = "https://fal.run/fal-ai/gemini-3-pro-image-preview"


def data_uri(path):
    mime = mimetypes.guess_type(path)[0] or "image/png"
    with open(path, "rb") as f:
        b = base64.b64encode(f.read()).decode()
    return f"data:{mime};base64,{b}"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--image", required=True)
    ap.add_argument("--prompt-file", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--aspect", default="21:9")
    ap.add_argument("--resolution", default="4K")
    ap.add_argument("--seed", type=int, default=None)
    a = ap.parse_args()

    key = os.environ.get("FAL_KEY")
    if not key:
        sys.exit("FAL_KEY not set")

    with open(a.prompt_file) as f:
        prompt = f.read().strip()

    payload = {
        "prompt": prompt,
        "image_urls": [data_uri(a.image)],
        "num_images": 1,
        "aspect_ratio": a.aspect,
        "resolution": a.resolution,
    }
    if a.seed is not None:
        payload["seed"] = a.seed

    req = urllib.request.Request(
        ENDPOINT,
        data=json.dumps(payload).encode(),
        headers={
            "Authorization": f"Key {key}",
            "Content-Type": "application/json",
            "Accept": "application/json",
        },
        method="POST",
    )
    try:
        with urllib.request.urlopen(req, timeout=300) as r:
            res = json.load(r)
    except urllib.error.HTTPError as e:
        body = e.read().decode(errors="replace")
        sys.exit(f"HTTP {e.code}: {body}")

    imgs = res.get("images") or []
    if not imgs:
        sys.exit(f"No images in response: {json.dumps(res)[:800]}")

    url = imgs[0]["url"]
    if url.startswith("data:"):
        b64 = url.split(",", 1)[1]
        raw = base64.b64decode(b64)
    else:
        with urllib.request.urlopen(url, timeout=300) as ir:
            raw = ir.read()
    with open(a.out, "wb") as f:
        f.write(raw)
    print(f"OK wrote {a.out} ({len(raw)} bytes); meta={json.dumps({k:v for k,v in res.items() if k!='images'})[:300]}")


if __name__ == "__main__":
    main()
