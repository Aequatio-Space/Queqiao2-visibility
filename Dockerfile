FROM docker.io/library/python:3.13-slim
LABEL "language"="python"
LABEL "framework"="flask"
WORKDIR /app
RUN apt-get update && apt-get install -y build-essential pkg-config clang && rm -rf /var/lib/apt/lists/*
# X11 开发库：pyvista/VTK 离屏渲染所需（容器内以 root 运行，无需 sudo）
RUN apt-get update && apt-get install -y libx11-dev libxrandr-dev libxinerama-dev libxcursor-dev libxi-dev && rm -rf /var/lib/apt/lists/*
RUN apt-get update && sudo apt-get install -y libegl1-mesa-dev libgles2-mesa-dev libglfw3-dev
RUN apt-get install -y libosmesa6-dev freeglut3-dev
COPY . .
RUN sed '/-e/d' requirements.txt | pip install -r /dev/stdin
RUN pip install -r requirements.txt
EXPOSE 8080
# 本仓库无 main.py，改用 gunicorn 启动（与 render.yaml/Procfile 一致；
# --workers 1 --threads 1 因 SPICE 非线程安全）
CMD ["/bin/bash", "-c", "_startup() { gunicorn viewer.server:app --bind 0.0.0.0:8080 --workers 1 --threads 1 --timeout 300; }; _startup"]