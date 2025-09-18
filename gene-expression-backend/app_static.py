from flask import Flask, render_template, jsonify
import os
import time
import fsspec
from pathlib import Path

from manifest_loader import get_manifest
from plotting_backend import get_interacting_partners

BASE_DIR = Path(__file__).resolve().parent

app = Flask(
	__name__,
	template_folder=str(BASE_DIR / "templates"),
	static_folder=str(BASE_DIR / "static"),
)

IMAGE_BASE_URL = os.getenv("IMAGE_BASE_URL", "")
USE_S3 = os.getenv("USE_S3", "false").lower() == "true"
S3_BUCKET = os.getenv("S3_BUCKET", "")

@app.route("/")
def home_static():
	# Provide a default S3 website style base if not set and in S3 mode
	base = IMAGE_BASE_URL
	if not base and USE_S3 and S3_BUCKET:
		base = f"https://{S3_BUCKET}.s3.amazonaws.com/data"
	return render_template("index_static.html", IMAGE_BASE_URL=base)

@app.route("/api/interacting-partners")
def api_partners():
	partners = get_interacting_partners()
	return {
		"interacting_partners": partners,
		"cytokines": sorted(partners.keys()),
		"total_pairs": sum(len(v) for v in partners.values()),
	}

@app.route("/api/tissues")
def api_tissues():
	# Return tissues without probing S3; assume available for static image flow
	items = {}
	manifest = get_manifest()
	for tissue, entry in manifest.items():
		key = entry.get("default")
		rec = {"key": key, "available": True}
		items[tissue] = rec
	return items

if __name__ == "__main__":
	print("🚀 Starting STATIC Flask server (serving pre-generated images)")
	print(f"   IMAGE_BASE_URL = {IMAGE_BASE_URL or '(derived)'}")
	app.run(debug=False, host="127.0.0.1", port=8002) 