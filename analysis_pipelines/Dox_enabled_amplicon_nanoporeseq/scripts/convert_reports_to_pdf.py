
import os
import glob
import subprocess

REPORT_DIR = "results/reports"
STYLE_FILE = "report_style.css"
MD_FILES = glob.glob(os.path.join(REPORT_DIR, "*.md"))

def convert_md_to_pdf(md_path):
    pdf_path = md_path.replace(".md", ".pdf")
    print(f"Converting {md_path} to {pdf_path}...")

    # Use pandoc with weasyprint and custom CSS
    cmd = [
        "pandoc",
        "-s", md_path,
        "-o", pdf_path,
        "--pdf-engine=weasyprint",
        "--css", STYLE_FILE,
        "--metadata", f"title={os.path.basename(md_path).replace('.md', '').replace('_', ' ').title()}"
    ]

    try:
        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode == 0:
            print(f"Successfully created {pdf_path}")
        else:
            print(f"Error converting {md_path}: {result.stderr}")
    except Exception as e:
        print(f"Failed to run pandoc for {md_path}: {e}")

def main():
    for md_file in MD_FILES:
        # Skip task.md if desired, but user said "all the md reports"
        convert_md_to_pdf(md_file)

if __name__ == "__main__":
    main()
