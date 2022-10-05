

def add_filter_args(parser):
    parser.add_argument('--min-hac', type=int, help="Min value for heavy atom count")
    parser.add_argument('--max-hac', type=int, help="Max value for heavy atom count")
    parser.add_argument('--min-rotb', type=int, help="Min value for rotatable bond count")
    parser.add_argument('--max-rotb', type=int, help="Max value for rotatable bond count")
    parser.add_argument('--min-rings', type=int, help="Min value for ring count")
    parser.add_argument('--max-rings', type=int, help="Max value for ring count")
    parser.add_argument('--min-aro-rings', type=int, help="Min value for aromatic ring count")
    parser.add_argument('--max-aro-rings', type=int, help="Max value for aromatic ring count")
    parser.add_argument('--min-chiral-centres', type=int, help="Min value for number of tetrahedral chiral centres")
    parser.add_argument('--max-chiral-centres', type=int, help="Max value for number of tetrahedral chiral centres")
    parser.add_argument('--min-undefined-chiral-centres', type=int,
                        help="Min value for number of undefined tetrahedral chiral centres")
    parser.add_argument('--max-undefined-chiral-centres', type=int,
                        help="Max value for number of undefined tetrahedral chiral centres")
    parser.add_argument('--min-sp3', type=int, help="Min value for SP3 count")
    parser.add_argument('--max-sp3', type=int, help="Max value for SP3 count")
    parser.add_argument('--min-logp', type=float, help="Min value for logP")
    parser.add_argument('--max-logp', type=float, help="Max value for logP")
    parser.add_argument('--min-tpsa', type=float, help="Min value for tpsa")
    parser.add_argument('--max-tpsa', type=float, help="Max value for tpsa")