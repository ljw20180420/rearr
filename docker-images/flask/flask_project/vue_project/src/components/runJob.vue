<script setup>
import { computed } from 'vue';
import axios from 'axios';
import { useNode, Handle, Position } from '@vue-flow/core';
import { validTargets } from './utils.js';

const base_url = import.meta.env.BASE_URL;
const { node, connectedEdges } = useNode();

node.active = computed(() => {
    for (let edge of connectedEdges.value) {
        if (edge.target === node.id && !edge.animated) {
            return false;
        }
    }
    return true;
});

node.data = computed(() => {
    let data = {};
    for (const edge of connectedEdges.value) {
        if (edge.source === node.id) {
            continue;
        }
        for (const name in edge.sourceNode.data) {
            data[name] = edge.sourceNode.data[name];
        }
    }
    return data;
})

function recursiveAssign(target, source) {
    for (const name in source) {
        if (typeof source[name] != 'object') {
            target[name] = source[name];
        } else if (name in target) {
            recursiveAssign(target[name], source[name]);
        }
    }
}

async function recursiveInspect(data) {
    if ('type' in data && typeof data.type == 'string') {
        if (data.taskId) {
            const response = await axios.get(base_url + "/inspect/" + data.taskId);
            if (response.status != 200) {
                alert(response.data);
                return false;
            }
            if (response.data != "SUCCESS") {
                alert(`job: ${param.taskId} is ${response.data}`);
                return false;
            }
        }
        return true;
    }
    for (const name in data) {
        if (!await recursiveInspect(data[name])) {
            return false;
        }
    }
    return true;
}

async function runJob() {
    try{
        if (!await recursiveInspect(node.data)) {
            return;
        }
        const response = await axios.put(`${base_url}/runJob/${node.id}`, node.data);
        if (response.status != 200) {
            alert(response.data);
            return;
        }
        for (let edge of connectedEdges.value) {
            if (edge.target === node.id) {
                continue;
            }
            recursiveAssign(edge.targetNode.data, response.data);
        }
    } catch(error) {
        alert(error);
    }
}
</script>

<template>
    <div :title="node.title" style="width: 200px; overflow: hidden;">
        <span>{{ node.id }} </span>
        <button @click="runJob" :disabled="!node.active">run</button>
        <Handle id="s" type="source" :position="Position.Right" :is-valid-connection="c => validTargets[c.source].includes(c.target)" />
        <Handle id="t" type="target" :position="Position.Left" :is-valid-connection="c => validTargets[c.source].includes(c.target)" />
    </div>
</template>