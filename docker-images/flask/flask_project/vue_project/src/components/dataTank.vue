<script setup>
import { ref, computed } from 'vue';
import axios from 'axios';
import { useNode, Handle, Position, useVueFlow, useHandleConnections } from '@vue-flow/core';
import { isValidConnection } from './utils.js'
const files = ref(null);
const percentage = ref(0);
const { node } = useNode();
const { findNode } = useVueFlow();

node.active = computed(() => {
    for (let obj of node.data) {
        if (obj.value == null || obj.value === "") {
            return false;
        } else if (Array.isArray(obj.value)) {
            for (let val of obj.value) {
                if (val == null) {
                    return false;
                }
            }
        }
    }
    return true;
})

useHandleConnections({
    type: "source",
    onDisconnect: (params) => {
        for (let param of params) {
            const target = findNode(param.target)
            for (let obj of node.data) {
                delete target.data.values[obj.name];
                if (Object.keys(obj).includes("taskId")) {
                    delete target.data.taskIds[obj.name];
                }
            }
        }
    }
});

async function upload(event, idx) {
    files.value = event.target.files;
    try{
        const response = await axios.putForm(base_url + "/upload", {
            files: files.value,
        },
        {
            onUploadProgress : progressEvent => {
                percentage.value = Math.round((progressEvent.loaded * 100) / progressEvent.total);
            }
        });
        node.data[idx].taskId = response.data.taskId;
        for (let i in node.data[idx].value) {
            node.data[idx].value[i] = response.data.value[i];
        }
    } catch(error) {
        alert(error);
    }
}

const base_url = import.meta.env.BASE_URL;
</script>

<template>
    <div :title="node.title">
        <div v-for="obj, i in node.data" :key="i">
            <span>{{ obj.name }}: </span>
            <a v-if="Object.keys(obj).includes('taskId')" :href="base_url + '/download/' + obj.taskId">{{ obj.taskId?.slice(0, 8) }}</a>
            <br>
            <input v-if="obj.type == 'value'" v-model="node.data[i].value">
            <div v-else-if="obj.type == 'file' || obj.type == 'files'" style="width: 200px; overflow: hidden;">
                <input :id="`${node.id}-input-file-${i}`" type="file" @change="upload($event, i)" :multiple="obj.type == 'files'">
                <div>{{ percentage }}</div>
            </div>
            <select v-else-if="obj.type == 'select'" v-model="node.data[i].value">
                <option v-for="option of node.data[i].options" :key="option">
                    {{ option }}
                </option>
            </select>
            <hr>
        </div>
        <Handle id="s" type="source" :position="Position.Right" :is-valid-connection="isValidConnection" />
        <Handle id="t" type="target" :position="Position.Left" :is-valid-connection="isValidConnection" />
    </div>
</template>