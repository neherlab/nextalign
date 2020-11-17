
dev:
	nodemon

dev-nowatch:
	CMAKE_BUILD_TYPE=Debug scripts/build_locally.sh


prod:
	CMAKE_BUILD_TYPE=Release scripts/build_locally.sh
