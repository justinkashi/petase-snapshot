#!/usr/bin/env python3
"""PETase paper search MCP server — no API key needed.

Self-contained RAG: pypdf for text extraction, sentence-transformers for
local embeddings (allenai-specter, science-domain), numpy for similarity search.
Claude Code synthesizes the final answer from returned passages.
"""

import pickle
import textwrap
from pathlib import Path

import numpy as np
from mcp.server.fastmcp import FastMCP

PAPER_DIR = Path(__file__).parent.parent / "petase_paperdb"
INDEX_PATH = Path(__file__).parent.parent / ".paperqa_index" / "index.pkl"
EMBED_MODEL = "allenai-specter"
CHUNK_SIZE = 600   # words per chunk
CHUNK_OVERLAP = 80

mcp = FastMCP("petase-papers")

# ── cached state ──────────────────────────────────────────────────────────────
_index: dict | None = None  # {"chunks": [...], "embeddings": np.ndarray, "sources": [...]}


def _get_model():
    from sentence_transformers import SentenceTransformer
    return SentenceTransformer(EMBED_MODEL)


def _chunk_text(text: str, source: str) -> list[dict]:
    words = text.split()
    chunks = []
    step = CHUNK_SIZE - CHUNK_OVERLAP
    for i in range(0, len(words), step):
        chunk_words = words[i : i + CHUNK_SIZE]
        if len(chunk_words) < 40:
            break
        chunks.append({"text": " ".join(chunk_words), "source": source})
    return chunks


def _extract_pdf(path: Path) -> str:
    from pypdf import PdfReader
    reader = PdfReader(str(path))
    return " ".join(
        page.extract_text() or "" for page in reader.pages
    )


def _load_index() -> dict | None:
    global _index
    if _index is not None:
        return _index
    if INDEX_PATH.exists():
        with open(INDEX_PATH, "rb") as f:
            _index = pickle.load(f)
    return _index


# ── tools ──────────────────────────────────────────────────────────────────────

@mcp.tool()
def build_index() -> str:
    """Build the PETase paper vector index (run once; ~5-15 min on CPU).

    Downloads allenai-specter model locally (~90 MB), processes all PDFs in
    petase_paperdb/, saves index to .paperqa_index/. No API key needed.
    """
    global _index

    model = _get_model()
    pdfs = sorted(PAPER_DIR.glob("**/*.pdf"))
    all_chunks, all_sources, added, failed = [], [], [], []

    for pdf in pdfs:
        try:
            text = _extract_pdf(pdf)
            chunks = _chunk_text(text, pdf.stem)
            all_chunks += [c["text"] for c in chunks]
            all_sources += [c["source"] for c in chunks]
            added.append(pdf.name)
        except Exception as e:
            failed.append(f"{pdf.name}: {e}")

    embeddings = model.encode(all_chunks, show_progress_bar=False, batch_size=32)

    _index = {"chunks": all_chunks, "embeddings": embeddings, "sources": all_sources}
    INDEX_PATH.parent.mkdir(parents=True, exist_ok=True)
    with open(INDEX_PATH, "wb") as f:
        pickle.dump(_index, f)

    msg = f"Indexed {len(added)} papers into {len(all_chunks)} chunks."
    if failed:
        msg += f"\nFailed: " + ", ".join(failed)
    return msg


@mcp.tool()
def query_papers(question: str, top_k: int = 8) -> str:
    """Search PETase literature and return relevant passages.

    Searches 46 papers using local vector similarity (no API key needed).
    Returns raw excerpts with source labels — Claude synthesizes the answer.
    Run build_index() first if not done yet.

    Args:
        question: Any biology/chemistry/engineering question about PETase
        top_k: Number of passages to return (default 8)
    """
    idx = _load_index()
    if idx is None:
        return "Index not built yet. Call build_index() first (takes ~5-15 min)."

    model = _get_model()
    q_emb = model.encode([question])
    scores = (idx["embeddings"] @ q_emb.T).squeeze()
    top_i = np.argsort(scores)[::-1][:top_k]

    parts = []
    for i in top_i:
        score = float(scores[i])
        source = idx["sources"][i]
        text = textwrap.fill(idx["chunks"][i], width=100)
        parts.append(f"[{source}] (score={score:.3f})\n{text}")

    return "\n\n---\n\n".join(parts)


@mcp.tool()
def list_papers() -> str:
    """List all papers in petase_paperdb/ and whether the index is built."""
    pdfs = sorted(PAPER_DIR.glob("**/*.pdf"))
    status = "index built" if INDEX_PATH.exists() else "index NOT built — run build_index()"
    lines = [f"{len(pdfs)} papers ({status}):"] + [f"  {p.name}" for p in pdfs]
    return "\n".join(lines)


if __name__ == "__main__":
    mcp.run(transport="stdio")
