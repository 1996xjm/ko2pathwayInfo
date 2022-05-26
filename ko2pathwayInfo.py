import json
import optparse
import pandas as pd
from io import StringIO

def transform(path):
    with open("./ko_index.json", "r") as f:
        json_dict = json.load(f)

    ko_info_strIO = StringIO()
    # with open("../data/kegg/ko_pathway_info.tsv", "w") as kpi_f:
    ko_info_strIO.write(f"user_id\tquery_ko\tpathway_count\tlevel1_pathway_id\tlevel1_pathway_name\tlevel2_pathway_id\tlevel2_pathway_name\tlevel3_pathway_id\tlevel3_pathway_name\tko\tko_name\tko_des\tec\n")
    with open(path, "r") as user_ko_f:
        for line in user_ko_f:
            split_res = line.split("\t")
            user_id = split_res[0].strip()
            user_ko = ""
            pathway_count = ""
            ko_info = ""
            if len(split_res) > 1:
                user_ko = split_res[1].strip()
                if user_ko in json_dict:
                    ko_info = json_dict[user_ko]["children"][-1]
                    pathway_count =  json_dict[user_ko]["count"]
                    for ko_info in json_dict[user_ko]["children"]:

                        ko_info_strIO.write(f"{user_id}\t{user_ko}\t{pathway_count}\t{ko_info}\n")
        # print(ko_info_strIO.getvalue())
        # 把读写位置恢复到起始，因为每次write后读写位置都在末尾，这样pandas读不到东西
        ko_info_strIO.seek(0)
        df = pd.read_table(ko_info_strIO)
        counts = df["level2_pathway_name"].value_counts()
        counts.name = "count"
        counts.index.name = "level2_pathway_name"
        res_json_dict = {
            "query_count":df.shape[0],
            "annotated_count":df.shape[0]-int(df["query_ko"].isnull().sum()),
            "level2_pathway_count_tsv":counts.to_csv(sep="\t"),
            "ko_info_tsv":ko_info_strIO.getvalue()
        }
        print(json.dumps(res_json_dict))


def transform2(path):
    """level3 id以逗号隔开"""
    with open("./ko_index.json", "r") as f:
        json_dict = json.load(f)

    ko_info_strIO = StringIO()
    # with open("../data/kegg/ko_pathway_info.tsv", "w") as kpi_f:
    ko_info_strIO.write(f"user_id\tquery_ko\tpathway_count\tlevel3_pathway_id\n")
    with open(path, "r") as user_ko_f:
        for line in user_ko_f:
            split_res = line.split("\t")
            user_id = split_res[0].strip()
            if len(split_res) > 1:
                user_ko = split_res[1].strip()
                if user_ko in json_dict:
                    pathway_count =  json_dict[user_ko]["count"]
                    level3_id = []
                    for ko_info in json_dict[user_ko]["children"]:
                        level3_id.append(ko_info.split("\t")[4])

                    ko_info_strIO.write(f"{user_id}\t{user_ko}\t{pathway_count}\t{','.join(level3_id)}\n")
        # print(ko_info_strIO.getvalue())
        # 把读写位置恢复到起始，因为每次write后读写位置都在末尾，这样pandas读不到东西
        ko_info_strIO.seek(0)
        df = pd.read_table(ko_info_strIO, index_col=0)

        df.to_csv("K_level3_id.tsv",sep="\t")




if __name__ == '__main__':
    usage = "python %prog -f/--filePath <ko list file>"
    parser = optparse.OptionParser(usage,description="Example: python %prog -f forward.fasta -r reverse.fasta -o assembly.fasta")  ## 写入上面定义的帮助信息
    parser.add_option('-f', '--filePath', dest='filePath', type='string', help='File path tsv file')
    options, args = parser.parse_args()
    transform(options.filePath)
    # transform2(options.filePath)
