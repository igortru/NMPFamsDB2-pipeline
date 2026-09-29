# NMPFamsDB2 pipeline tools, versions as tested (linux/amd64; blast 2.17.0 has no linux-aarch64 build)
# env "nmpfams": search, clustering, alignment, structure search
# env "colabfold": structure prediction (CPU jaxlib from conda; GPU needs a CUDA jaxlib)
FROM --platform=linux/amd64 mambaorg/micromamba:1.5.10

RUN micromamba create -y -n nmpfams -c conda-forge -c bioconda \
        python=3.11 \
        pandas=3.0.6 \
        numpy=2.4.6 \
        scikit-learn=1.9.1 \
        hmmer=3.4 \
        tantan=51 \
        last=1654 \
        diamond=2.2.8 \
        mmseqs2=18.8cc5c \
        mafft=7.526 \
        hhsuite=3.3.0 \
        foldseek=10.941cd33 \
        blast=2.17.0 \
        lbzip2=2.5 \
    && micromamba create -y -n colabfold -c conda-forge -c bioconda \
        python=3.11 \
        colabfold=1.5.5 \
    && micromamba clean -a -y

COPY --chown=$MAMBA_USER:$MAMBA_USER aligner.py trimmer.py redundancy_removal.py parser.sh commands.md README.md /opt/NMPFamsDB2/

ENV PATH=/opt/conda/envs/nmpfams/bin:$PATH \
    MPLCONFIGDIR=/tmp/mpl
WORKDIR /data

# colabfold: micromamba run -n colabfold colabfold_batch ...
