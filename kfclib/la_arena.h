#ifndef LA_ARENA_H
#define LA_ARENA_H

struct la_arena
{
	size_t size;
	void *start;
	void *end;
	void *cur;

};

struct la_arena *la_arena_create(size_t nbytes);
void *la_arena_alloc(struct la_arena *arena, size_t nbytes);
void la_arena_free(struct la_arena *arena);
void la_arena_reset(struct la_arena *arena);

#endif
