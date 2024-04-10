/* Copyright 2024 Informatics Matters Ltd.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
*/

params.digits = 6
params.prefix = 'x'
params.suffix = '.txt'
params.header = false // if true then skip the first line
params.chunk_size = 1000

process split_txt {

    container 'informaticsmatters/vs-prep:stable'

    input:
    file inputs

    output:
    file params.prefix + '*' + params.suffix

    """
    ${params.header ? 'tail -n +2' : 'cat'} $inputs |
    split -l $params.chunk_size -d -a $params.digits --additional-suffix $params.suffix - $params.prefix
    """
}