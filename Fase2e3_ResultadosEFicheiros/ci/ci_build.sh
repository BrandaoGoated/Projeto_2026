PREFIX="/opt/conda"
SRC="/app"

mkdir -p "${PREFIX}/bin"

cp "${SRC}/m_party.py" "${PREFIX}/bin/"
cp -r "${SRC}/workflow" "${PREFIX}/bin/"
cp -r "${SRC}/resources" "${PREFIX}/bin/"
cp -r "${SRC}/config" "${PREFIX}/bin/"
cp "${SRC}/ci/sequences_PET.fasta" "${PREFIX}/bin/"

chmod +x "${PREFIX}/bin/m_party.py"

# alias
ln -s /opt/conda/bin/m_party.py /opt/conda/bin/m-party