<script setup>
import { computed } from 'vue';
import { useNode, Handle, Position } from '@vue-flow/core';
import { validTargets } from './utils.js';
import recursiveParams from './recursiveParams.vue';

function recursiveCheckActive(data) {
    if ('type' in data && typeof data.type == 'string') {
        return data.value != null;
    }
    for (const name in data) {
        if (!recursiveCheckActive(data[name])) {
            return false;
        }
    }
    return true;
}

const { node } = useNode();

node.active = computed(() => recursiveCheckActive(node.data));
</script>

<template>
    <div :title="node.title" style="width: 200px; overflow: hidden;">
        <div style="display: flex; flex-direction: row; align-items: center; justify-content: center;">
            {{ node.id }}
        </div>
        <hr width="100%">
        <recursiveParams v-model="node.data" />
        <Handle id="s" type="source" :position="Position.Right" :is-valid-connection="c => validTargets[c.source].includes(c.target)" />
        <Handle id="t" type="target" :position="Position.Left" :is-valid-connection="c => validTargets[c.source].includes(c.target)" />
    </div>
</template>