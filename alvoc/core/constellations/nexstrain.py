# nexstrain.py

from typing import Iterator, Tuple, List, Dict, Any
from alvoc.core.constellations.core import DataSource

class NextstrainSource(DataSource):
    """
    Nextstrain-specific data source for JSON tree dumps.
    """

    def fetch(self, source: str) -> Dict[str, Any]:
        import requests
        headers = {"Accept": "application/vnd.nextstrain.dataset.main+json"}
        resp = requests.get(source, headers=headers)
        resp.raise_for_status()
        return resp.json()

    def records(self, raw_data: Dict[str, Any]) -> Iterator[Tuple[str, List[str]]]:
        from collections import deque

        # 1. find the root of the Nextstrain tree
        root = raw_data["tree"]

        # 2. first pass: collect every mutation for bitmask indexing
        muts = set()
        queue = deque([root])
        while queue:
            node = queue.popleft()
            muts.update(
                node.get("branch_attrs", {})
                    .get("mutations", {})
                    .get("nuc", [])
            )
            queue.extend(node.get("children", []))

        mutations_list = sorted(muts)
        idx_map = {m: i for i, m in enumerate(mutations_list)}

        # 3. second pass: traverse tree carrying an inherited bitmask
        queue = deque([(root, 0)])
        while queue:
            node, inh_mask = queue.popleft()

            # build mask for this branch’s mutations
            curr_mask = 0
            for m in node.get("branch_attrs", {}).get("mutations", {}).get("nuc", []):
                curr_mask |= 1 << idx_map[m]

            # union with inherited mask
            all_mask = inh_mask | curr_mask

            # drop true reversions
            tmp = curr_mask
            while tmp:
                lowbit = tmp & -tmp
                i = lowbit.bit_length() - 1
                mut = mutations_list[i]
                ref, pos, alt = mut[0], mut[1:-1], mut[-1]
                rev = f"{alt}{pos}{ref}"
                if rev in idx_map and ((inh_mask >> idx_map[rev]) & 1):
                    all_mask &= ~(1 << i)
                tmp ^= lowbit

            # extract clade name and collect this node’s mutation list
            clade = node.get("node_attrs", {}) \
                        .get("clade_membership", {}) \
                        .get("value")
            if clade:
                muts_here: List[str] = []
                mask = all_mask
                while mask:
                    lowbit = mask & -mask
                    idx = lowbit.bit_length() - 1
                    muts_here.append(mutations_list[idx])
                    mask ^= lowbit

                # **yield** each record for calculate_profiles()
                yield clade, muts_here

            # enqueue children with the updated mask
            for child in node.get("children", []):
                queue.append((child, all_mask))
