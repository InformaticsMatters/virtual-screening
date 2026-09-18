# Copyright 2025 Informatics Matters Ltd.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Mapping a molecule digest onto the sharded directory layout used by the
enumeration/conformer pipeline.

This lives here rather than in dm_job_utilities: the shared package
deliberately dropped get_path_from_digest as obsolete, but the layout it
describes is still what prepare_enum_conf_lists.py writes and
assemble_conformers.py reads back.
"""

default_num_chars = 2
default_num_levels = 2


def get_path_from_digest(
    digest, num_chars=default_num_chars, num_levels=default_num_levels
):
    """Split a digest into the directory components that locate its files,
    e.g. 'abcdef...' -> ['ab', 'cd'] for the default 2 chars over 2 levels.
    """
    parts = []
    start = 0
    for _ in range(0, num_levels):
        end = start + num_chars
        parts.append(digest[start:end])
        start = start + num_chars
    return parts
