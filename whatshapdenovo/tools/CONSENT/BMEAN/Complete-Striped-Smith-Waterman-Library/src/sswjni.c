#include <jni.h>
#include <stdio.h>
#include "ssw.h"

jstring s_align_cigar_to_jstring(JNIEnv* env, s_align* align) {
	jstring jstrBuf = NULL;
	if (align->cigar != NULL && align->cigarLen > 0) {
		char* buffer = (char *)malloc(16 * align->cigarLen + 1);
		buffer[0] = '\0';
		char* currentBufferPosition = buffer;
		int i;
		for (i = 0; i < align->cigarLen; i++) {
			int charsPrinted = sprintf(currentBufferPosition, "%d%c", cigar_int_to_len(align->cigar[i]), cigar_int_to_op(align->cigar[i]));
			currentBufferPosition += charsPrinted;
		}
		jstrBuf = (*env)->NewStringUTF(env, buffer);
		free(buffer);
	}
	return jstrBuf;
}
jobject s_align_to_ssw_Alignment(JNIEnv* env, s_align* align) {
	if (align == NULL) return NULL;
	jclass clazz = (*env)->FindClass(env, "ssw/Alignment");
	jmethodID constructor = (*env)->GetMethodID(env, clazz, "<init>", "(SSIIIIILjava/lang/String;)V");
	jobject result = (*env)->NewObject(env, clazz, constructor,
		align->score1,
		align->score2,
		align->ref_begin1,
		align->ref_end1,
		align->read_begin1,
		align->read_end1,
		align->ref_end2,
		s_align_cigar_to_jstring(env, align));
	return result;
}
JNIEXPORT jobject JNICALL Java_ssw_Aligner_align(JNIEnv* env, jclass cls,
		jbyteArray read, jbyteArray matrix, jint matrixSize, jint score_size,
		jbyteArray ref,
		jint gapOpen,
		jint gapExtend,
		jint flag,
		jshort filters,
		jint filterd,
		jint maskLen) {
	jbyte* readPtr = (*env)->GetByteArrayElements(env, read, NULL);
	jsize readLen = (*env)->GetArrayLength(env, read);
	jbyte* matrixPtr = (*env)->GetByteArrayElements(env, matrix, NULL);
	/*jsize matrixLen = (*env)->GetArrayLength(env, matrix);*/
	jbyte* refPtr = (*env)->GetByteArrayElements(env, ref, NULL);
	jsize refLen = (*env)->GetArrayLength(env, ref);
	s_profile* profile = ssw_init(readPtr, readLen, matrixPtr, matrixSize, (int8_t)score_size);
	s_align* align = ssw_align(profile, refPtr, refLen, gapOpen, gapExtend, flag, filters, filterd, maskLen);
	jobject jalignment = s_align_to_ssw_Alignment(env, align);
	align_destroy(align);
	init_destroy(profile);
	(*env)->ReleaseByteArrayElements(env, read, readPtr, JNI_ABORT);
	(*env)->ReleaseByteArrayElements(env, matrix, matrixPtr, JNI_ABORT);
	(*env)->ReleaseByteArrayElements(env, ref, refPtr, JNI_ABORT);
	return jalignment;
}

