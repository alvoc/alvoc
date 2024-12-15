import aiohttp
import aiofiles
import json
from pathlib import Path
import asyncio


async def download_lineage_list(input_url: str, outdir: Path):
    async with aiohttp.ClientSession() as session:
        async with session.get(input_url) as response:
            if response.status != 200:
                raise ConnectionError("Failed to download lineage list.")
            text = await response.text()

    lines = text.splitlines()

    output_file = outdir / "lineages.txt"
    async with aiofiles.open(output_file, "w") as file:
        await file.write("\n".join(lines))
    values = [
        line.split()[0]
        for line in lines[1:]
        if line.strip() and not line.startswith("*")
    ]

    return values  # Directly return the list of lineages


async def fetch_lineage_data(session, lineage, base_url, start_date):
    url = f"{base_url}/sample/aggregated?dateFrom={start_date}&pangoLineage={lineage}"
    async with session.get(url) as response:
        if response.status != 200:
            print(f"Error fetching data for lineage {lineage}")
            return None
        data = await response.json()

    if data["data"][0]["count"] > 100:
        mutations_url = f"{base_url}/sample/nucleotideMutations?pangoLineage={lineage}"
        async with session.get(mutations_url) as mutations_response:
            if mutations_response.status == 200:
                mutations_data = await mutations_response.json()
                mutations = [
                    m["mutation"]
                    for m in mutations_data["data"]
                    if m["proportion"] > 0.90 and not m["mutation"].endswith("-")
                ]
                return {
                    "lineage": lineage,
                    "label": f"{lineage}-like",
                    "description": f"{lineage} lineage defining mutations",
                    "sources": [],
                    "tags": [lineage],
                    "sites": mutations,
                    "note": "Unique mutations for sublineage",
                }
            else:
                print(f"Error fetching mutations for lineage {lineage}")
    return None


async def make_constellations(input_url, outdir, base_url):
    outdir.mkdir(parents=True, exist_ok=True)
    lineages = await download_lineage_list(input_url, outdir)
    start_date = "2021-01-01"
    all_lineage_data = {}

    async with aiohttp.ClientSession() as session:
        tasks = [
            fetch_lineage_data(session, lineage, base_url, start_date)
            for lineage in lineages
        ]
        results = await asyncio.gather(*tasks)

        for result in results:
            if result:
                all_lineage_data[result["lineage"]] = result

    # Writing all data into a single JSON file
    output_file = outdir / "covid_constellations.json"
    async with aiofiles.open(output_file, "w") as file:
        await file.write(json.dumps(all_lineage_data, indent=4))


# Run the asynchronous function
async def main():
    outdir = Path("constellations")
    input_url = "https://raw.githubusercontent.com/cov-lineages/pango-designation/master/lineage_notes.txt"
    base_url = "https://lapis.cov-spectrum.org/open/v2"
    await make_constellations(input_url, outdir, base_url)


asyncio.run(main())
