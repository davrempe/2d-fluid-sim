#version 440

in vec2 vPos;
in vec3 vColor;

uniform mat4 MVP;

out vec3 fColor;

void main() {
	fColor = vColor;
	// add MVP
	gl_Position = MVP * vec4(vPos.x, vPos.y, 0.0, 1.0);
}