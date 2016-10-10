#include "LoadProgram.h"
#include <stdio.h>
#include <stdlib.h>
#include <gl/glut.h>

static char* readShaderSource(const char* shaderFile)
{
	FILE* fp = fopen(shaderFile, "r");

	if (fp == NULL) { return NULL; }

	fseek(fp, 0L, SEEK_END);
	long size = ftell(fp);

	fseek(fp, 0L, SEEK_SET);
	char* buf = new char[size + 1];
	fread(buf, 1, size, fp);

	buf[size] = '\0';
	fclose(fp);

	return buf;
}

GLuint LoadProgram(const char* vShaderFile, const char* fShaderFile)
{
	// read in shader contents
	const char* vShader = readShaderSource(vShaderFile);
	const char* fShader = readShaderSource(fShaderFile);

	GLuint shader, program;
	GLint completed;

	program = glCreateProgram();

	// load and compile vertex shader
	if (vShader != NULL) {
		shader = glCreateShader(GL_VERTEX_SHADER);
		glShaderSource(shader, 1, &vShader, NULL);
		glCompileShader(shader);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &completed);

		if (!completed) {
			GLint len;
			char* msg;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
			msg = (char*)malloc(len);
			glGetShaderInfoLog(shader, len, &len, msg);
			fprintf(stderr, "Vertex shader compilation failure:\n%s\n", msg);
			free(msg);
			glDeleteProgram(program);

			exit(EXIT_FAILURE);
		}

		glAttachShader(program, shader);
		
	}

	// load and compile fragment shader
	if (fShader != NULL) {
		shader = glCreateShader(GL_FRAGMENT_SHADER);
		glShaderSource(shader, 1, &fShader, NULL);
		glCompileShader(shader);
		glGetShaderiv(shader, GL_COMPILE_STATUS, &completed);

		if (!completed) {
			GLint len;
			char* msg;
			glGetShaderiv(shader, GL_INFO_LOG_LENGTH, &len);
			msg = (char*)malloc(len);
			glGetShaderInfoLog(shader, len, &len, msg);
			fprintf(stderr, "Fragment shader compilation failure:\n%s\n", msg);
			free(msg);
			glDeleteProgram(program);
			exit(EXIT_FAILURE);
		}

		glAttachShader(program, shader);
	}

	// link program
	glLinkProgram(program);
	glGetProgramiv(program, GL_LINK_STATUS, &completed);

	if (!completed) {
		GLint len;
		char* msg;
		glGetProgramiv(program, GL_INFO_LOG_LENGTH, &len);
		msg = (char*)malloc(len);
		glGetProgramInfoLog(program, len, &len, msg);
		fprintf(stderr, "Program link failure:\n%s\n", msg);
		free(msg);
		glDeleteProgram(program);
		exit(EXIT_FAILURE);
	}

	return program;

}