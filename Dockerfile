# Reproducible environment: Ubuntu 24.04 + Python deps (pinned) + FEBio 4.13 from source with MKL PARDISO.
FROM ubuntu:24.04
ENV DEBIAN_FRONTEND=noninteractive PIP_BREAK_SYSTEM_PACKAGES=1
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential cmake git ca-certificates python3 python3-pip \
      libglu1-mesa libxrender1 libxcursor1 libxft2 libxinerama1 libgl1 libosmesa6 xvfb \
      texlive-latex-base texlive-latex-recommended texlive-latex-extra && rm -rf /var/lib/apt/lists/*
WORKDIR /opt/rlfea
COPY scripts/build_febio.sh scripts/build_febio.sh
RUN bash scripts/build_febio.sh && ln -s /opt/febio-4.13/bin/febio4 /usr/local/bin/febio4
COPY requirements.lock.txt .
RUN python3 -m pip install -r requirements.lock.txt
COPY . .
RUN python3 -m pip install -e . --no-deps
CMD ["rlfea", "--help"]
