FROM ubuntu:24.04 AS build

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
    build-essential \
    ca-certificates \
    cmake \
    ninja-build \
    pkg-config \
    libgsl-dev \
    python3 \
    python3-pip \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /opt/hyperiso
COPY . .

RUN cmake -S Hyperiso/Hyperiso/core -B build \
      -G Ninja \
      -DCMAKE_BUILD_TYPE=Release \
      -DBUILD_WITH_CLI=ON \
      -DENABLE_TESTS=OFF \
 && cmake --build build --parallel 2 \
 && cmake --install build --prefix /usr/local

CMD ["hyperiso-ui", "--help"]
