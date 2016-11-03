
in vec2 vPos;
in vec3 vColor;

uniform mat4 MVP;

out vec3 fColor;

void main() {
	fColor = vColor;

	gl_Position = MVP * vec4(vPos, 0.0, 1.0);
}