#!/usr/bin/env python3
import os
import argparse
from pathlib import Path
import time
import sys

# Ensure project root (gene-expression-backend) is on sys.path
CURRENT_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = CURRENT_DIR.parent
if str(PROJECT_ROOT) not in sys.path:
	sys.path.insert(0, str(PROJECT_ROOT))

import scanpy as sc

from manifest_loader import get_manifest, dataset_key_for
from data_access import load_dataset_for
from plotting_backend import plot_histogram_for_pair

try:
	import boto3
	HAVE_BOTO3 = True
except Exception:
	HAVE_BOTO3 = False

from plotting_backend import get_interacting_partners


def ensure_dir(p: Path):
	p.mkdir(parents=True, exist_ok=True)


def save_fig(fig, path: Path):
	ensure_dir(path.parent)
	fig.savefig(str(path), format="png", bbox_inches="tight", dpi=100)
	print(f"   wrote {path}")
	# Close the figure to free memory during long runs
	import matplotlib.pyplot as plt
	plt.close(fig)


def generate_for_tissue(tissue: str, out_root: Path):
	print(f"🧩 Generating for {tissue}")
	adata = load_dataset_for(tissue)
	adata.obs['tissue_type'] = tissue
	pairs = []
	partners = get_interacting_partners()
	for g1, recs in partners.items():
		for g2 in recs:
			pairs.append((g1, g2))

	tissue_dir = out_root / tissue
	for (g1, g2) in pairs:
		fig = plot_histogram_for_pair(adata, tissue, g1, g2)
		path = tissue_dir / f"{g1}_{g2}.png"
		save_fig(fig, path)


def upload_directory(local_dir: Path, bucket: str, prefix: str = "data/"):
	if not HAVE_BOTO3:
		raise SystemExit("boto3 not installed; pip install boto3")
	s3 = boto3.client('s3')
	for path in local_dir.rglob('*.png'):
		key = f"{prefix}{path.relative_to(local_dir).as_posix()}"
		s3.upload_file(str(path), bucket, key, ExtraArgs={'ContentType': 'image/png', 'ACL': 'public-read'})
		print(f"   uploaded s3://{bucket}/{key}")


def main():
	ap = argparse.ArgumentParser()
	ap.add_argument('--out', default='static_plots', help='output directory for PNGs')
	ap.add_argument('--tissue', action='append', help='limit to specific tissue(s)')
	ap.add_argument('--upload', action='store_true', help='upload to S3 (requires AWS creds)')
	ap.add_argument('--bucket', default=os.getenv('S3_BUCKET', ''), help='S3 bucket for upload')
	ap.add_argument('--prefix', default='data/', help='S3 key prefix for uploads')
	args = ap.parse_args()

	out_root = Path(args.out).resolve()
	ensure_dir(out_root)
	print(f"Output directory: {out_root}")

	manifest = get_manifest()
	tissues = args.tissue or list(manifest.keys())
	for t in tissues:
		generate_for_tissue(t, out_root)

	if args.upload:
		if not args.bucket:
			raise SystemExit('--bucket is required for --upload')
		upload_directory(out_root, args.bucket, args.prefix)

if __name__ == '__main__':
	main() 