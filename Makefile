dev:
	@$(MAKE) --no-print-directory dev-impl

dev-impl:
	@nodemon

dev-nowatch:
	CMAKE_BUILD_TYPE=Debug scripts/build_locally.sh

prod:
	CMAKE_BUILD_TYPE=Release scripts/build_locally.sh

format:
	@scripts/format.sh
