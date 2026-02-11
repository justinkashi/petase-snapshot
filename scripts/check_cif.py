# Usage mode 1 (QC): python scripts/check_cif.py --cif file1.cif file2.cif
# Usage mode 2 (Align): python scripts/check_cif.py --cif file.cif --fasta seq.fasta --chain-id A
import argparse
from Bio import SeqIO, pairwise2
from Bio.PDB import MMCIFParser, PDBParser, PPBuilder
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from Bio.SeqUtils import seq1

parser = argparse.ArgumentParser(description="Compare one or more mmCIF/PDB sequences to a FASTA sequence.")
parser.add_argument("--cif", required=True, nargs="+", help="Path(s) to input .cif/.mmcif/.pdb file(s)")
parser.add_argument("--fasta", help="Path to input .fasta file (optional)")
parser.add_argument("--chain-id", default="", help="Chain ID (e.g., A) to compare only that chain")
args = parser.parse_args()

cif_paths = args.cif
fasta_path = args.fasta
chain_id = args.chain_id.strip()

fasta_seq = None
if fasta_path:
    fasta_seq = str(next(SeqIO.parse(fasta_path, "fasta")).seq).upper()

def align_stats(a, b):
    aln = pairwise2.align.globalms(a, b, 2, -1, -10, -0.5, one_alignment_only=True)[0]
    matches = sum(1 for x, y in zip(aln.seqA, aln.seqB) if x == y and x != '-' and y != '-')
    aligned = sum(1 for x, y in zip(aln.seqA, aln.seqB) if x != '-' and y != '-')
    pid = matches / aligned if aligned else 0
    cov = aligned / len(a) if len(a) else 0
    return pid, cov


def is_cif(path):
    p = path.lower()
    return p.endswith(".cif") or p.endswith(".mmcif")


def load_structure(path):
    parser = MMCIFParser(QUIET=True) if is_cif(path) else PDBParser(QUIET=True)
    return parser.get_structure("s", path)


def canonical_cif_seqs(path):
    d = MMCIF2Dict(path)
    if "_entity_poly_seq.entity_id" not in d:
        return None
    eids = d["_entity_poly_seq.entity_id"]
    mons = d["_entity_poly_seq.mon_id"]
    tmp = {}
    for eid, mon in zip(eids, mons):
        tmp.setdefault(eid, []).append(mon)
    seqs = {}
    for eid, monlist in tmp.items():
        seqs[f"entity_{eid}"] = seq1(" ".join(monlist), custom_map={"MSE": "M"})
    return seqs


def chain_seqs(structure):
    ppb = PPBuilder()
    seqs = {}
    for chain in structure[0]:
        seq = "".join(str(pp.get_sequence()) for pp in ppb.build_peptides(chain))
        if seq:
            seqs[f"chain_{chain.id}"] = seq
    return seqs

for cif_path in cif_paths:
    print(f"\nCIF: {cif_path}")

    if fasta_seq:
        print("FASTA length:", len(fasta_seq))
        seqs = {}
        if chain_id:
            # Force using the observed sequence from a specific chain
            structure = load_structure(cif_path)
            target = None
            for chain in structure[0]:
                if chain.id == chain_id:
                    target = chain
                    break
            if target is None:
                raise ValueError(f"Chain '{chain_id}' not found in {cif_path}")
            seq = "".join(str(pp.get_sequence()) for pp in PPBuilder().build_peptides(target))
            if seq:
                seqs[f"chain_{target.id}"] = seq
        else:
            # Prefer full canonical sequences if present (mmCIF only)
            if is_cif(cif_path):
                seqs = canonical_cif_seqs(cif_path) or {}
            if not seqs:
                # Fallback: observed chain sequences
                structure = load_structure(cif_path)
                seqs = chain_seqs(structure)

        for name, seq in seqs.items():
            pid, cov = align_stats(fasta_seq, seq)
            print(name, "len", len(seq), "identity", f"{pid*100:.1f}%", "coverage", f"{cov*100:.1f}%")
    else:
        # QC-only mode
        structure = load_structure(cif_path)
        models = len(structure)
        model = structure[0]
        print("models:", models)
        print("chains:", len(model))
        warn = []
        if models > 1:
            warn.append("multiple models")
        for chain in model:
            residues = list(chain)
            aa_res = [r for r in residues if r.id[0] == " "]
            het_res = [r for r in residues if r.id[0] != " "]
            atom_count = sum(len(r) for r in residues)
            resnums = [r.id[1] for r in aa_res if isinstance(r.id[1], int)]
            min_res = min(resnums) if resnums else None
            max_res = max(resnums) if resnums else None
            if len(chain.id) != 1:
                warn.append(f"chain id '{chain.id}' not 1 char")
            print(
                f"chain {chain.id} residues {len(residues)} aa {len(aa_res)} "
                f"hetero {len(het_res)} atoms {atom_count} "
                f"resnum_range {min_res}..{max_res}"
            )
        if warn:
            print("warnings:", ", ".join(sorted(set(warn))))
