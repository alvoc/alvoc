

<div style="display: flex; justify-content: center; align-items: center; text-align: center; flex-direction: column;">
  <div style="display: flex; align-items: center;">
    <img src="https://raw.githubusercontent.com/alvoc/alvoc/main/docs/assets/icon.svg" alt="Logo" width="200" height="auto">
    <span style="font-size: 100px; color: #4e8ada;"> Alvoc </span>
  </div>

  <p><em>Abundance learning for variants of concern</em></p>

  <div style="text-align: center;">
    <p>
      <a href="https://github.com/alvoc/alvoc/issues" target="_blank">Report Bug</a>
      |
      <a href="https://github.com/alvoc/alvoc/issues" target="_blank">Request Feature</a>
    </p>
    <p>
      <a href="https://alvoc.github.io" target="_blank"> Documentation </a>
    </p>
  </div> 
</div>

---

## Overview

Alvoc is a tool for estimating the abundace of variants of concern from wastewater sequencing data. Its predecessor is [Alcov](https://github.com/Ellmen/alcov), an abundance learning tool for SARS-CoV-2 variants.

In addition alvoc can also be used for:

* Determining the frequency of mutations of interest in BAM files
* Converting nucleotide and amino acid mutations for a particular virus
* Comparing amplicon GC content with its read depth (as a measure of degredation)

## Installation

### With Pip

Pip is the default package installer for Python.

1. Install [Python](https://www.python.org/downloads/). We recommend using [pyenv](https://github.com/pyenv/pyenv) for python version management.

2. Install the latest version from pypi.

```console
pip install alvoc
```

### With UV

UV is a fast, all-in-one Python tool for dependency, version, and project management, replacing pip, poetry, and more.

1. Install [UV](https://github.com/astral-sh/uv).
2. Install the latest version from pypi.

```console
uv add alvoc
```

Yes, including a Docker image under the installation section makes sense, especially for users who prefer containerized environments or need to integrate your tool into CI/CD pipelines. Here's how you can include it:

---

### With Docker

For users who prefer containerized environments, an official Docker image is available on GitHub Container Registry (GHCR).

1. Ensure you have [Docker](https://www.docker.com) installed and running.
2. Pull the latest Docker image:

```console
docker pull ghcr.io/alvoc/alvoc:latest
```

3. Run the container:

```console
docker run --rm -it ghcr.io/alvoc/alvoc:latest
```

---

## What is Alvoc?  

Alvoc is a command-line tool designed to analyze viral sequencing data with a focus on identifying **variants of concern** and their defining mutations. With the ongoing rise of infectious diseases like **SARS-CoV-2**, monitoring and understanding viral evolution has become critical for researchers, public health officials, and laboratories.  

### Why Use Alvoc?  

1. **Lineage Identification**  
    - Alvoc helps determine which viral lineages are present in your sequencing data. For example, it can identify emerging lineages like **B.1.1.7 (Alpha)** or **P.1 (Gamma)** for SARS-CoV-2.  
    - Use Alvoc when you need to:  
        - Compare samples to known lineage-defining mutations.  
        - Detect specific lineages in samples from an outbreak.  
        - Include or exclude specific lineages (whitelisting/blacklisting).  

2. **Mutation Detection**  
    - Alvoc allows you to detect mutations in sequencing data with customizable depth thresholds.  
    - Use this feature to:  
        - Find known or novel mutations in viral samples.  
        - Focus on specific mutations like **S:N501Y** or **E484K** that are critical for vaccine efficacy or transmissibility.  

3. **Amplicon Quality Assessment**  
    - Viral sequencing often relies on **amplicon-based approaches**. Alvoc provides tools to assess amplicon metrics like **coverage** and **GC content**.  
    - This helps ensure high-quality sequencing data for downstream analyses.  

### When to Use Alvoc?  

You should use Alvoc if you are:  

- A **researcher** analyzing viral sequencing data to study genetic variation and evolution.  
- A **public health analyst** monitoring variants of concern during outbreaks.  
- A **bioinformatician** processing large-scale sequencing data to identify mutations and lineages.  
- A **laboratory scientist** validating amplicon-based sequencing data for accuracy.  

Alvoc simplifies these workflows into a set of streamlined commands, saving time and improving consistency across analyses.  

---


# License
This project is licensed under the terms of the MIT license.