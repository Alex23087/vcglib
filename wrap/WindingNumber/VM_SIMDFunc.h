#pragma once
#ifndef __SSE__
#ifndef __VM_SIMDFunc__
#define __VM_SIMDFunc__

#include "SYS_Types.h"

#include <cmath>

namespace igl { namespace FastWindingNumber {

struct v4si {
	int32 v[4];
};

struct v4sf {
	float v[4];
};

static SYS_FORCE_INLINE v4sf V4SF(const v4si &v) {
	static_assert(sizeof(v4si) == sizeof(v4sf) && alignof(v4si) == alignof(v4sf), "v4si and v4sf must be compatible");
	return *(const v4sf*)&v;
}

static SYS_FORCE_INLINE v4si V4SI(const v4sf &v) {
	static_assert(sizeof(v4si) == sizeof(v4sf) && alignof(v4si) == alignof(v4sf), "v4si and v4sf must be compatible");
	return *(const v4si*)&v;
}

static SYS_FORCE_INLINE int32 conditionMask(bool c) {
	return c ? int32(0xFFFFFFFF) : 0;
}

static SYS_FORCE_INLINE v4sf
VM_SPLATS(float f) {
	return v4sf{{f, f, f, f}};
}

static SYS_FORCE_INLINE v4si
VM_SPLATS(uint32 i) {
	return v4si{{int32(i), int32(i), int32(i), int32(i)}};
}

static SYS_FORCE_INLINE v4si
VM_SPLATS(int32 i) {
	return v4si{{i, i, i, i}};
}

static SYS_FORCE_INLINE v4sf
VM_SPLATS(float a, float b, float c, float d) {
	return v4sf{{a, b, c, d}};
}

static SYS_FORCE_INLINE v4si
VM_SPLATS(uint32 a, uint32 b, uint32 c, uint32 d) {
	return v4si{{int32(a), int32(b), int32(c), int32(d)}};
}

static SYS_FORCE_INLINE v4si
VM_SPLATS(int32 a, int32 b, int32 c, int32 d) {
	return v4si{{a, b, c, d}};
}

static SYS_FORCE_INLINE v4si
VM_LOAD(const int32 v[4]) {
	return v4si{{v[0], v[1], v[2], v[3]}};
}

static SYS_FORCE_INLINE v4sf
VM_LOAD(const float v[4]) {
	return v4sf{{v[0], v[1], v[2], v[3]}};
}


static inline v4si VM_ICMPEQ(v4si a, v4si b) {
	return v4si{{
		conditionMask(a.v[0] == b.v[0]),
		conditionMask(a.v[1] == b.v[1]),
		conditionMask(a.v[2] == b.v[2]),
		conditionMask(a.v[3] == b.v[3])
	}};
}

static inline v4si VM_ICMPGT(v4si a, v4si b) {
	return v4si{{
		conditionMask(a.v[0] > b.v[0]),
		conditionMask(a.v[1] > b.v[1]),
		conditionMask(a.v[2] > b.v[2]),
		conditionMask(a.v[3] > b.v[3])
	}};
}

static inline v4si VM_ICMPLT(v4si a, v4si b) {
	return v4si{{
		conditionMask(a.v[0] < b.v[0]),
		conditionMask(a.v[1] < b.v[1]),
		conditionMask(a.v[2] < b.v[2]),
		conditionMask(a.v[3] < b.v[3])
	}};
}

static inline v4si VM_IADD(v4si a, v4si b) {
	return v4si{{
		(a.v[0] + b.v[0]),
		(a.v[1] + b.v[1]),
		(a.v[2] + b.v[2]),
		(a.v[3] + b.v[3])
	}};
}

static inline v4si VM_ISUB(v4si a, v4si b) {
	return v4si{{
		(a.v[0] - b.v[0]),
		(a.v[1] - b.v[1]),
		(a.v[2] - b.v[2]),
		(a.v[3] - b.v[3])
	}};
}

static inline v4si VM_OR(v4si a, v4si b) {
	return v4si{{
		(a.v[0] | b.v[0]),
		(a.v[1] | b.v[1]),
		(a.v[2] | b.v[2]),
		(a.v[3] | b.v[3])
	}};
}

static inline v4si VM_AND(v4si a, v4si b) {
	return v4si{{
		(a.v[0] & b.v[0]),
		(a.v[1] & b.v[1]),
		(a.v[2] & b.v[2]),
		(a.v[3] & b.v[3])
	}};
}

static inline v4si VM_ANDNOT(v4si a, v4si b) {
	return v4si{{
		((~a.v[0]) & b.v[0]),
		((~a.v[1]) & b.v[1]),
		((~a.v[2]) & b.v[2]),
		((~a.v[3]) & b.v[3])
	}};
}

static inline v4si VM_XOR(v4si a, v4si b) {
	return v4si{{
		(a.v[0] ^ b.v[0]),
		(a.v[1] ^ b.v[1]),
		(a.v[2] ^ b.v[2]),
		(a.v[3] ^ b.v[3])
	}};
}

static SYS_FORCE_INLINE int
VM_EXTRACT(const v4si v, int index) {
	return v.v[index];
}

static SYS_FORCE_INLINE float
VM_EXTRACT(const v4sf v, int index) {
	return v.v[index];
}

static SYS_FORCE_INLINE v4si
VM_INSERT(v4si v, int32 value, int index) {
	v.v[index] = value;
	return v;
}

static SYS_FORCE_INLINE v4sf
VM_INSERT(v4sf v, float value, int index) {
	v.v[index] = value;
	return v;
}

static inline v4si VM_CMPEQ(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] == b.v[0]),
		conditionMask(a.v[1] == b.v[1]),
		conditionMask(a.v[2] == b.v[2]),
		conditionMask(a.v[3] == b.v[3])
	}};
}

static inline v4si VM_CMPNE(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] != b.v[0]),
		conditionMask(a.v[1] != b.v[1]),
		conditionMask(a.v[2] != b.v[2]),
		conditionMask(a.v[3] != b.v[3])
	}};
}

static inline v4si VM_CMPGT(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] > b.v[0]),
		conditionMask(a.v[1] > b.v[1]),
		conditionMask(a.v[2] > b.v[2]),
		conditionMask(a.v[3] > b.v[3])
	}};
}

static inline v4si VM_CMPLT(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] < b.v[0]),
		conditionMask(a.v[1] < b.v[1]),
		conditionMask(a.v[2] < b.v[2]),
		conditionMask(a.v[3] < b.v[3])
	}};
}

static inline v4si VM_CMPGE(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] >= b.v[0]),
		conditionMask(a.v[1] >= b.v[1]),
		conditionMask(a.v[2] >= b.v[2]),
		conditionMask(a.v[3] >= b.v[3])
	}};
}

static inline v4si VM_CMPLE(v4sf a, v4sf b) {
	return v4si{{
		conditionMask(a.v[0] <= b.v[0]),
		conditionMask(a.v[1] <= b.v[1]),
		conditionMask(a.v[2] <= b.v[2]),
		conditionMask(a.v[3] <= b.v[3])
	}};
}

static inline v4sf VM_ADD(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] + b.v[0]),
		(a.v[1] + b.v[1]),
		(a.v[2] + b.v[2]),
		(a.v[3] + b.v[3])
	}};
}

static inline v4sf VM_SUB(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] - b.v[0]),
		(a.v[1] - b.v[1]),
		(a.v[2] - b.v[2]),
		(a.v[3] - b.v[3])
	}};
}

static inline v4sf VM_NEG(v4sf a) {
	return v4sf{{
		(-a.v[0]),
		(-a.v[1]),
		(-a.v[2]),
		(-a.v[3])
	}};
}

static inline v4sf VM_MUL(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] * b.v[0]),
		(a.v[1] * b.v[1]),
		(a.v[2] * b.v[2]),
		(a.v[3] * b.v[3])
	}};
}

static inline v4sf VM_DIV(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] / b.v[0]),
		(a.v[1] / b.v[1]),
		(a.v[2] / b.v[2]),
		(a.v[3] / b.v[3])
	}};
}

static inline v4sf VM_MADD(v4sf a, v4sf b, v4sf c) {
	return v4sf{{
		(a.v[0] * b.v[0]) + c.v[0],
		(a.v[1] * b.v[1]) + c.v[1],
		(a.v[2] * b.v[2]) + c.v[2],
		(a.v[3] * b.v[3]) + c.v[3]
	}};
}

static inline v4sf VM_ABS(v4sf a) {
	return v4sf{{
		(a.v[0] < 0) ? -a.v[0] : a.v[0],
		(a.v[1] < 0) ? -a.v[1] : a.v[1],
		(a.v[2] < 0) ? -a.v[2] : a.v[2],
		(a.v[3] < 0) ? -a.v[3] : a.v[3]
	}};
}

static inline v4sf VM_MAX(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] < b.v[0]) ? b.v[0] : a.v[0],
		(a.v[1] < b.v[1]) ? b.v[1] : a.v[1],
		(a.v[2] < b.v[2]) ? b.v[2] : a.v[2],
		(a.v[3] < b.v[3]) ? b.v[3] : a.v[3]
	}};
}

static inline v4sf VM_MIN(v4sf a, v4sf b) {
	return v4sf{{
		(a.v[0] > b.v[0]) ? b.v[0] : a.v[0],
		(a.v[1] > b.v[1]) ? b.v[1] : a.v[1],
		(a.v[2] > b.v[2]) ? b.v[2] : a.v[2],
		(a.v[3] > b.v[3]) ? b.v[3] : a.v[3]
	}};
}

static inline v4sf VM_INVERT(v4sf a) {
	return v4sf{{
		(1.0f/a.v[0]),
		(1.0f/a.v[1]),
		(1.0f/a.v[2]),
		(1.0f/a.v[3])
	}};
}

static inline v4sf VM_SQRT(v4sf a) {
	return v4sf{{
		std::sqrt(a.v[0]),
		std::sqrt(a.v[1]),
		std::sqrt(a.v[2]),
		std::sqrt(a.v[3])
	}};
}

static inline v4si VM_INT(v4sf a) {
	return v4si{{
		int32(a.v[0]),
		int32(a.v[1]),
		int32(a.v[2]),
		int32(a.v[3])
	}};
}

static inline v4sf VM_IFLOAT(v4si a) {
	return v4sf{{
		float(a.v[0]),
		float(a.v[1]),
		float(a.v[2]),
		float(a.v[3])
	}};
}

static SYS_FORCE_INLINE void VM_P_FLOOR() {}

static SYS_FORCE_INLINE int32 singleIntFloor(float f) {
	// Casting to int32 usually truncates toward zero, instead of rounding down,
	// so subtract one if the result is above f.
	int32 i = int32(f);
	i -= (float(i) > f);
	return i;
}
static inline v4si VM_FLOOR(v4sf a) {
	return v4si{{
		singleIntFloor(a.v[0]),
		singleIntFloor(a.v[1]),
		singleIntFloor(a.v[2]),
		singleIntFloor(a.v[3])
	}};
}

static SYS_FORCE_INLINE void VM_E_FLOOR() {}

static SYS_FORCE_INLINE bool vm_allbits(v4si a) {
	return (
		(a.v[0] == -1) && 
		(a.v[1] == -1) && 
		(a.v[2] == -1) && 
		(a.v[3] == -1)
	);
}

int SYS_FORCE_INLINE _mm_movemask_ps(const v4si& v) {
	return (
		int(v.v[0] < 0) |
		(int(v.v[1] < 0)<<1) |
		(int(v.v[2] < 0)<<2) |
		(int(v.v[3] < 0)<<3)
	);
}

int SYS_FORCE_INLINE _mm_movemask_ps(const v4sf& v) {
	// Use std::signbit just in case it needs to distinguish between +0 and -0
	// or between positive and negative NaN values (e.g. these could really
	// be integers instead of floats).
	return (
		int(std::signbit(v.v[0])) |
		(int(std::signbit(v.v[1]))<<1) |
		(int(std::signbit(v.v[2]))<<2) |
		(int(std::signbit(v.v[3]))<<3)
	);
}
}}
#endif
#endif
