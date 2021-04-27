# Computer graphics and interaction (DH2323) labs

This repository is the implementation of the rasterization track of DH2323 computer graphics and interaction course at KTH.

All labs use SDL2 instead of SDL. Implementation of SDL2 wrapper was received from: https://github.com/lemonad/DH2323-Skeleton

### Build
```
cd Lab-*
mkdir build && cd build
cmake ..
make
```

### Lab 1 - Interpolation
#### Color
Illustrates interpolation of colors onto the screen

##### Run
Within build directory, find _Color_
```
./Color
```

#### Starfield
Illustrates a classic star-wars like star field with perspective projection
##### Run
Within build directory, find _Starfield_
```
./Starfield
```

### Lab 2 - Raytracing
Illustrates a simple ray-traced Cornell box with 5 walls and 2 blocks including diffuse and specular lighting and shadows.

##### Run
Within build directory, find _Raytracing_
```
./Raytracing
```

### Lab 3 - Rasterization
Illustrates a simple rasterized Cornell box with 5 walls and 2 blocks including depth testing and diffuse and specular lighting.

##### Run
Within build directory, find _Rasterization_
```
./Rasterization
```
