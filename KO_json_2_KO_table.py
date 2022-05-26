import json
import re
import pandas as pd

def transformJson2table():
    with open("./ko00001.json") as f:
        ko_map_data = json.load(f)

    with open("KEGG_pathway_ko.tsv", "w") as oh:
        line = "level1_pathway_id\tlevel1_pathway_name\tlevel2_pathway_id\tlevel2_pathway_name"
        line += "\tlevel3_pathway_id\tlevel3_pathway_name\tko\tko_name\tko_des\tec\n"
        oh.write(line)
        level1_list = ko_map_data["children"]
        level1_list.reverse()
        for level1 in level1_list:
            m = re.match(r"(\S+)\s+([\S\w\s]+)", level1["name"])
            level1_pathway_id = m.groups()[0].strip()
            level1_pathway_name = m.groups()[1].strip()
            for level2 in level1["children"]:
                m = re.match(r"(\S+)\s+([\S\w\s]+)", level2["name"])
                level2_pathway_id = m.groups()[0].strip()
                level2_pathway_name = m.groups()[1].strip()
                for level3 in level2["children"]:
                    m = re.match(r"(\S+)\s+([^\[]*)", level3["name"])
                    level3_pathway_id = m.groups()[0].strip()
                    level3_pathway_name = m.groups()[1].strip()
                    if "children" in level3:
                        for ko in level3["children"]:
                            m = re.match(r"(\S+)\s+(.+?);\s+(.*)( \[EC:.*\])", ko["name"])
                            if m is not None:
                                ko_id = m.groups()[0].strip()
                                ko_name = m.groups()[1].strip()
                                ko_des = m.groups()[2].strip()
                                ec = m.groups()[3]
                            else:
                                m = re.match(r"(\S+)\s+(.+?);\s+(.*)", ko["name"])
                                if m is not None:
                                    ko_id = m.groups()[0].strip()
                                    ko_name = m.groups()[1].strip()
                                    ko_des = m.groups()[2].strip()
                                    ec = "-"
                                else:
                                    print(ko["name"])
                            line = level1_pathway_id + "\t" + level1_pathway_name + "\t" + level2_pathway_id + "\t" + level2_pathway_name
                            line += "\t" + level3_pathway_id + "\t" + level3_pathway_name + "\t" + ko_id + "\t" + ko_name + "\t" + ko_des + "\t" + ec + "\n"
                            oh.write(line)

                    else:
                        print("no child",level3_pathway_name)


def koTable2JsonIndex():
    """ko表转为以ko为key的json"""
    df = pd.read_table("./KEGG_pathway_ko.tsv", dtype="string")
    json_dict = {}
    for key,group_df in df.groupby("ko"):
        group_df = group_df.set_index("level1_pathway_id")
        info_list = group_df.to_csv(sep="\t").split("\n")[1:-1]
        json_dict[key] = {
            "children":info_list,
            "count":len(info_list)
        }
    with open("./ko_index.json", "w") as f:
        json.dump(json_dict, f, indent=4)


if __name__ == '__main__':
    koTable2JsonIndex()

