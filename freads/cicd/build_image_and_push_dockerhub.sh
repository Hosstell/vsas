set -e

docker build -t freads .

echo "Deploying latest version"
docker image tag freads vsas/freads:latest
docker image push vsas/freads:latest

VERSION=$(cat version.txt)
echo "Deploying ${VERSION} version"
docker image tag freads vsas/freads:${VERSION}
docker image push vsas/freads:${VERSION}
