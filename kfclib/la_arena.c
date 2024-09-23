#include <stdint.h>
#include <stdlib.h>
#include "la_arena.h"

#define ARCH_MAX_ALIGN (sizeof(void *))

static void *mem_align(struct la_arena *arena);


struct la_arena *la_arena_create(size_t nbytes)
{
	struct la_arena *res = malloc(sizeof(*res));
	if(res == NULL)
		return NULL;

	res->start = malloc(nbytes);
	if(res->start == NULL) {
		free(res);
		return NULL;
	}

	uintptr_t end = (uintptr_t)res->start + nbytes;
	res->end = (void *)end; 
	res->cur = res->start;

	return res;
}

void la_arena_reset(struct la_arena *arena)
{
	arena->cur = arena->start;
}

void la_arena_free(struct la_arena *arena)
{
	if(arena == NULL)
		return;
	
	free(arena->start);
	free(arena);
}


void *la_arena_alloc(struct la_arena *arena, size_t nbytes)
{
	void *start = mem_align(arena);	
	if(start == NULL)
		return NULL;

	uintptr_t end = nbytes + (uintptr_t)start;
	if(end >= (uintptr_t)arena->end)
		return NULL;	
	
	arena->cur = (void *)end; 
	return start;	
}


static void *mem_align(struct la_arena *arena)
{
	uintptr_t cur = (uintptr_t)arena->cur;
	uintptr_t offset = cur % ARCH_MAX_ALIGN;
	cur += offset;

	if(cur >= (uintptr_t)arena->end)
		return NULL;
	
	return (void *)cur;
}


