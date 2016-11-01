
in vec2 vPos;
in vec3 vColor;

uniform mat4 mProj;

out vec3 fColor;

void main() {
	fColor = vColor;

	gl_Position = mProj * vec4(vPos, 0.0, 1.0);
}