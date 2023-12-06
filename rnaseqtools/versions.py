def get_versions():
    return versions[0]["number"]


versions = [
    {
        "number": "0.0.3",
        "features": [
            "1. Fix some spelling mistakes",
            "2. Finish the readme of Count2TMM and GetGeneLength",
            "3. Fit the new version of ToolBox",
        ],
    },
    {
        "number": "0.0.2",
        "features": [
            "1. Add blat mapping for DenovoCount2RefCount",
            "2. Rearrange the code",
            "3. Finish the readme of DenovoCount2RefCount",
        ],
    },
    {
        "number": "0.0.1",
        "features": [
            "1. Separate code about RNA-seq from ToolBiox",
        ],
    },
]