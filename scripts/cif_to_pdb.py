import argparse
from pathlib import Path

import gemmi


def chain_name(chain):
    return getattr(chain, "name", getattr(chain, "id", ""))


def main():
    parser = argparse.ArgumentParser(description="Convert mmCIF/.cif to PDB with optional chain filtering.")
    parser.add_argument("--cif", required=True, help="Path to input .cif/.mmcif file")
    parser.add_argument("--out", default="", help="Path to output .pdb (default: <input>.pdb)")
    parser.add_argument("--chain-id", default="", help="Keep only this chain ID (e.g., A)")
    parser.add_argument(
        "--protein-only",
        action="store_true",
        help="Remove ligands and waters (keep only polymer residues)",
    )
    parser.add_argument("--remove-waters", action="store_true", help="Remove waters")
    parser.add_argument(
        "--remove-ligands",
        action="store_true",
        help="Remove ligands (also removes waters in gemmi)",
    )
    args = parser.parse_args()

    cif_path = Path(args.cif)
    out_path = Path(args.out) if args.out else cif_path.with_suffix(".pdb")

    st = gemmi.read_structure(str(cif_path))
    model = st[0]

    if args.chain_id:
        for ch in list(model):
            if chain_name(ch) != args.chain_id:
                model.remove_chain(chain_name(ch))
        if len(model) == 0:
            raise ValueError(f"Chain '{args.chain_id}' not found in {cif_path}")

    if args.protein_only:
        st.remove_ligands_and_waters()
    else:
        if args.remove_waters:
            st.remove_waters()
        if args.remove_ligands:
            st.remove_ligands_and_waters()

    st.write_pdb(str(out_path))
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
