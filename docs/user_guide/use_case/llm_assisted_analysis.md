# Using LLM Engines with scvi-tools

Large language models (LLMs) can significantly lower the barrier to using scvi-tools by helping researchers write code, choose models, tune parameters, and troubleshoot analyses through natural language. This page covers how to leverage five popular AI platforms—Claude, ChatGPT, OpenClaw, Gemini, and BioMNI—to get the most out of scvi-tools.

---

## Claude (Anthropic)

Claude is a general-purpose AI assistant from Anthropic. For scvi-tools users, Claude offers a dedicated **scvi-tools Skill Bundle**—a curated set of skills covering the full scvi-tools ecosystem.

### scvi-tools Skill Bundle

The skill bundle gives Claude deep knowledge of scvi-tools workflows and includes guidance for:

- **Batch integration**: scVI and scArches
- **Cell type annotation**: SCANVI and CellAssign
- **Spatial analysis**: DestVI, Tangram, Cell2location, Stereoscope
- **Epigenetic data**: PeakVI and scBasset
- **Multimodal integration**: TotalVI (CITE-seq) and MultiVI (RNA+ATAC)
- **Perturbation studies**: contrastiveVI

Each skill covers recommended workflows, parameter guidance, and troubleshooting tips.

### Installation

**Claude Code users:**
```bash
/plugin install scvi-tools@life-sciences
```

**Claude.ai users:**
Organization admins can upload the skill ZIP via *Admin Settings > Skills*. Individual users can upload via *Settings > Capabilities > Skills*. Download instructions and the skill ZIP are provided in the [Anthropic tutorial](https://claude.com/resources/tutorials/how-to-use-the-scvi-tools-bioinformatics-skill-bundle-with-claude).

Once installed, you can ask questions like:
> "I have 10x Chromium data from 3 donors with different sequencing depths. Which scvi-tools model should I use for integration, and what batch key should I set?"

See the full tutorial at [Anthropic's scvi-tools Skill Bundle guide](https://claude.com/resources/tutorials/how-to-use-the-scvi-tools-bioinformatics-skill-bundle-with-claude).

---

## ChatGPT (OpenAI)

ChatGPT can assist with scvi-tools through two complementary routes: custom GPTs and MCP (Model Context Protocol) tool integrations.

### Custom GPTs

OpenAI's GPT Store hosts community-built GPTs specialized in single-cell analysis. For example, the [Scanpy – Your Single-Cell RNA-seq Data Analyst](https://chatgpt.com/g/g-GKNExWk2P-scanpy-your-single-cell-rna-seq-data-analyst) GPT is configured to assist with scanpy-based workflows, which pair naturally with scvi-tools preprocessing pipelines.

You can use such GPTs to:
- Walk through an end-to-end scRNA-seq analysis
- Get scvi-tools code snippets for common operations
- Debug errors from scvi-tools model training

### MCP (Model Context Protocol) Tool Use

OpenAI supports [tool use via the API](https://developers.openai.com/api/docs/guides/tools/), which enables agents to call external functions—including Python code execution. This makes it possible to build automated pipelines where ChatGPT generates and runs scvi-tools code on your data.

A simple agent prompt example:
> "Load the AnnData file at `data/pbmc.h5ad`, run scVI with 2 batches defined by `adata.obs['batch']`, train for 400 epochs, and return the UMAP coordinates."

With tool use enabled, ChatGPT can generate the code and invoke a Python execution environment to produce results.

---

## OpenClaw

[OpenClaw](https://lobehub.com/skills/k-dense-ai-claude-scientific-skills-scvi-tools) (available via the LobeHub market) provides an installable skill focused on scvi-tools for use with Claude-based agents. It is optimized for researchers who need rigorous statistical frameworks and multi-batch integration.

### Installation

```bash
# Register your agent (one-time)
npx -y @lobehub/market-cli register --name "YourName" --source open-claw

# Install the scvi-tools skill
npx -y @lobehub/market-cli skills install k-dense-ai-claude-scientific-skills-scvi-tools
```

After installation, read the `SKILL.md` file in the extracted directory for usage instructions. The skill covers:

- Probabilistic batch correction and dataset alignment
- Multi-modal analysis (CITE-seq, spatial, multiome)
- Uncertainty quantification in differential expression
- Cell annotation with transfer learning

This skill is best suited for users who want a lightweight Claude-compatible skill without requiring the full Claude.ai platform.

---

## Gemini (Google)

Gemini is Google's general-purpose LLM, accessible via [Google AI Studio](https://aistudio.google.com) and the Gemini API. While there is no dedicated scvi-tools skill for Gemini, it is effective for assisted code generation, debugging, and conceptual guidance when working with scvi-tools.

### General Use

Gemini can help with scvi-tools through natural language prompting for:
- Generating scvi-tools setup and training code
- Explaining model outputs and hyperparameters
- Suggesting appropriate models for your data type

**Example prompt in AI Studio:**
> "Write Python code to train an scVI model on an AnnData object with a batch column called `sample_id`, then extract the latent embedding and run a UMAP."

You can paste scvi-tools error messages, documentation excerpts, or code snippets directly into the chat to get targeted assistance.

---

## BioMNI (Stanford)

[BioMNI](https://biomni.stanford.edu/) is a general-purpose biomedical AI agent from Stanford, designed to autonomously execute research tasks across diverse biomedical subfields. It has a native integration with the scverse ecosystem—including scvi-tools—announced in 2025.

### Integration with scverse and scvi-tools

BioMNI understands biological context and can orchestrate multi-step pipelines across scverse packages (Scanpy, scvi-tools, Squidpy, Pertpy) from plain-language instructions. Crucially, all agent-generated code is packaged as reproducible Jupyter notebooks.

**Example prompts:**
> "Run QC and normalization, integrate my three batches using scVI, cluster the cells, and annotate cell types using SCANVI."

> "Cluster cells and identify marker genes for each cluster."

BioMNI handles parameter selection, dependency management, and returns documented, reproducible results—no manual coding required.

### Access

- **Web platform**: [biomni.stanford.edu](https://biomni.stanford.edu/) — interactive, no setup required
- **Open-source**: [github.com/snap-stanford/Biomni](https://github.com/snap-stanford/Biomni) — self-hosted deployment

BioMNI is particularly well-suited for biologists who want to run complete single-cell and spatial workflows without writing code, while still producing reproducible, shareable analyses.

---

## Summary

| Platform | Best For | scvi-tools Integration |
|---|---|---|
| **Claude** | Guided workflows, parameter tuning, troubleshooting | Dedicated skill bundle with full model coverage |
| **ChatGPT** | Code generation, custom GPTs, agentic pipelines | Custom GPTs + MCP tool use |
| **OpenClaw** | Lightweight Claude-based skill, CLI install | Installable scvi-tools skill via LobeHub |
| **Gemini** | General code assistance, AI Studio prompting | General LLM assistance; no dedicated skill |
| **BioMNI** | End-to-end automated scverse pipelines | Native scverse/scvi-tools integration |

Each platform offers a different trade-off between ease of use, customization, and depth of scvi-tools knowledge. For users who primarily want guidance and code examples, Claude's skill bundle or BioMNI provide the deepest integration. For users building custom pipelines or agentic workflows, ChatGPT's MCP tool use or BioMNI's open-source deployment offer the most flexibility.

:::{note}
LLM-generated code should always be reviewed before running on important data. Check that model parameters, batch keys, and data shapes match your specific dataset.
:::
