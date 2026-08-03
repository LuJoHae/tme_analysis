import argparse
from pathlib import Path
from PIL import Image

def combine_horizontally(img1_path: Path, img2_path: Path, out_path: Path):
    try:
        img1 = Image.open(img1_path)
        img2 = Image.open(img2_path)
    except Exception as e:
        print(f"Error opening images: {e}")
        return

    # Resize img2 to match img1 height
    if img1.height != img2.height:
        new_width = int(img2.width * (img1.height / img2.height))
        img2 = img2.resize((new_width, img1.height), Image.Resampling.LANCZOS)
        
    combined = Image.new('RGB', (img1.width + img2.width, img1.height), (255, 255, 255))
    combined.paste(img1, (0, 0))
    combined.paste(img2, (img1.width, 0))
    
    out_path.parent.mkdir(parents=True, exist_ok=True)
    combined.save(out_path)
    print(f"Saved combined image to {out_path}")

def main():
    parser = argparse.ArgumentParser(description="Combine plots")
    parser.add_argument("--img1", required=True)
    parser.add_argument("--img2", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()
    
    combine_horizontally(Path(args.img1), Path(args.img2), Path(args.out))

if __name__ == "__main__":
    main()
