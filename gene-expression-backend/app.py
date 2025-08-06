# app.py
from flask import Flask, request, send_file
from plotting_backend import load_adata, plot_histogram_for_pair
import io
import os
import subprocess
import matplotlib.pyplot as plt

H5AD_FILENAME = "final_object.h5ad"
GDRIVE_FILE_ID = "11Un1mYmyvNWNlAnu2vW3ERMD0NdY7qJs"

# Download and unzip the file if it doesn't exist
if not os.path.exists(H5AD_FILENAME):
    print("📦 Downloading final_object.zip from Google Drive...")
    subprocess.run([
        "gdown", "--id", GDRIVE_FILE_ID, "--output", "final_object.zip"
    ], check=True)

    print("📂 Extracting final_object.h5ad from zip...")
    subprocess.run(["unzip", "-o", "final_object.zip"])

# Load the AnnData file once at startup
print("🧬 Loading AnnData...")
adata = load_adata(H5AD_FILENAME)

app = Flask(__name__)

@app.route("/")
def home():
    return "🎉 Flask server is running!"

@app.route("/generate-plot", methods=["GET"])
def generate_plot():
    tissue = request.args.get("tissue")
    gene1 = request.args.get("gene1")
    gene2 = request.args.get("gene2")

    if not all([tissue, gene1, gene2]):
        return "Missing required parameters", 400

    try:
        fig = plot_histogram_for_pair(adata, tissue, gene1, gene2)

        buf = io.BytesIO()
        fig.savefig(buf, format="png", bbox_inches="tight")
        buf.seek(0)
        plt.close(fig)

        return send_file(buf, mimetype="image/png")

    except Exception as e:
        return f"Error: {str(e)}", 500
